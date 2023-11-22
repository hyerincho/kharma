/* 
 *  File: imex_step.cpp
 *  
 *  BSD 3-Clause License
 *  
 *  Copyright (c) 2020, AFD Group at UIUC
 *  All rights reserved.
 *  
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *  
 *  1. Redistributions of source code must retain the above copyright notice, this
 *     list of conditions and the following disclaimer.
 *  
 *  2. Redistributions in binary form must reproduce the above copyright notice,
 *     this list of conditions and the following disclaimer in the documentation
 *     and/or other materials provided with the distribution.
 *  
 *  3. Neither the name of the copyright holder nor the names of its
 *     contributors may be used to endorse or promote products derived from
 *     this software without specific prior written permission.
 *  
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 *  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 *  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 *  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 *  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 *  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#include "kharma_driver.hpp"

#include "decs.hpp"

//Packages
#include "b_cd.hpp"
#include "b_cleanup.hpp"
#include "b_ct.hpp"
#include "b_flux_ct.hpp"
#include "electrons.hpp"
#include "grmhd.hpp"
#include "wind.hpp"
// Other headers
#include "boundaries.hpp"
#include "flux.hpp"
#include "kharma.hpp"
#include "resize_restart.hpp"
#include "implicit.hpp"

#include <parthenon/parthenon.hpp>
#include <interface/update.hpp>
#include <amr_criteria/refinement_package.hpp>

TaskCollection KHARMADriver::MakeImExTaskCollection(BlockList_t &blocks, int stage)
{
    // Reminder that this list is created BEFORE any of the list contents are run!
    // Prints or function calls here will likely not do what you want: instead, add to the list by calling tl.AddTask()

    TaskCollection tc;
    TaskID t_none(0);

    // Which packages we've loaded affects which tasks we'll add to the list
    auto& pkgs         = blocks[0]->packages.AllPackages();
    auto& driver_pkg   = pkgs.at("Driver")->AllParams();
    const bool use_electrons = pkgs.count("Electrons");
    const bool use_b_cleanup = pkgs.count("B_Cleanup");
    const bool use_b_ct = pkgs.count("B_CT");
    const bool use_implicit = pkgs.count("Implicit");
    const bool use_jcon = pkgs.count("Current");
    const bool use_linesearch = (use_implicit) ? pkgs.at("Implicit")->Param<bool>("linesearch") : false;

    // Allocate/copy the things we need
    // TODO these can now be reduced by including the var lists/flags which actually need to be allocated
    // TODO except the Copy they can be run on step 1 only
    if (stage == 1) {
        auto &base = pmesh->mesh_data.Get();
        // Fluxes
        pmesh->mesh_data.Add("dUdt");
        for (int i = 1; i < integrator->nstages; i++)
            pmesh->mesh_data.Add(integrator->stage_name[i]);
        // Preserve state for time derivatives if we need to output current
        if (use_jcon) {
            pmesh->mesh_data.Add("preserve");
            // Above only copies on allocate -- ensure we copy every step
            Copy<MeshData<Real>>({Metadata::Cell}, base.get(), pmesh->mesh_data.Get("preserve").get());
        }
        if (use_implicit) {
            // When solving, we need a temporary copy with any explicit updates,
            // but not overwriting the beginning- or mid-step values
            pmesh->mesh_data.Add("solver");
            if (use_linesearch) {
                // Need an additional state for linesearch
                pmesh->mesh_data.Add("linesearch");
            }
        }
    }

    static std::vector<std::string> sync_vars;
    if (sync_vars.size() == 0) {
        // Build the universe of variables to let Parthenon see when exchanging boundaries.
        // This is built to exclude incidental variables like B field initialization stuff, EMFs, etc.
        using FC = Metadata::FlagCollection;
        auto sync_flags = FC({Metadata::GetUserFlag("Primitive"), Metadata::Conserved, Metadata::Face}, true);
        sync_vars = KHARMA::GetVariableNames(&(pmesh->packages), sync_flags);
    }
    // We'll only ever sync the current stage "final"
    //pmesh->mesh_data.AddShallow("sync", integrator->stage_name[stage], sync_vars);

    // Big synchronous region: get & apply fluxes to advance the fluid state
    // num_partitions is nearly always 1
    const int num_partitions = pmesh->DefaultNumPartitions();
    TaskRegion &single_tasklist_per_pack_region = tc.AddRegion(num_partitions);
    for (int i = 0; i < num_partitions; i++) {
        auto &tl = single_tasklist_per_pack_region[i];
        // Container names: 
        // '_full_step_init' refers to the fluid state at the start of the full time step (Si in iharm3d)
        // '_sub_step_init' refers to the fluid state at the start of the sub step (Ss in iharm3d)
        // '_sub_step_final' refers to the fluid state at the end of the sub step (Sf in iharm3d)
        // '_flux_src' refers to the mesh object corresponding to -divF + S
        // '_solver' refers to the fluid state passed to the Implicit solver. At the end of the solve
        // '_linesearch' refers to the fluid state updated while performing a linesearch in the solver
        // copy P and U from solver state to sub_step_final state.
        auto &md_full_step_init = pmesh->mesh_data.GetOrAdd("base", i);
        auto &md_sub_step_init  = pmesh->mesh_data.GetOrAdd(integrator->stage_name[stage - 1], i);
        auto &md_sub_step_final = pmesh->mesh_data.GetOrAdd(integrator->stage_name[stage], i);
        auto &md_flux_src       = pmesh->mesh_data.GetOrAdd("dUdt", i);
        // Normally we put explicit update in md_solver, then add implicitly-evolved variables and copy back.
        // If we're not doing an implicit solve at all, just write straight to sub_step_final
        std::shared_ptr<MeshData<Real>> &md_solver = (use_implicit) ? pmesh->mesh_data.GetOrAdd("solver", i) : md_sub_step_final;
        auto &md_sync = pmesh->mesh_data.AddShallow("sync"+integrator->stage_name[stage]+std::to_string(i), md_sub_step_final, sync_vars);

        // Start receiving flux corrections and ghost cells
        auto t_start_recv_bound = tl.AddTask(t_none, parthenon::StartReceiveBoundBufs<parthenon::BoundaryType::any>, md_sub_step_final);
        auto t_start_recv_flux = t_start_recv_bound;
        if (pmesh->multilevel || use_b_ct)
            t_start_recv_flux = tl.AddTask(t_none, parthenon::StartReceiveFluxCorrections, md_sub_step_init);
        
        // Calculate the flux of each variable through each face
        // This reconstructs the primitives (P) at faces and uses them to calculate fluxes
        // of the conserved variables (U) through each face.
        const KReconstruction::Type& recon = driver_pkg.Get<KReconstruction::Type>("recon");
        auto t_fluxes = KHARMADriver::AddFluxCalculations(t_start_recv_bound, tl, recon, md_sub_step_init.get());

        // Any package modifications to the fluxes.  e.g.:
        // 1. CT calculations for B field transport
        // 2. Zero fluxes through poles
        // etc 
        auto t_fix_flux = tl.AddTask(t_fluxes, Packages::FixFlux, md_sub_step_init.get());

        // If we're in AMR, correct fluxes from neighbors
        auto t_flux_bounds = t_fix_flux;
        if (pmesh->multilevel || use_b_ct) {
            auto t_emf = t_flux_bounds;
            if (use_b_ct) {
                // Pull out a container of only EMF to synchronize
                auto &md_emf_only = pmesh->mesh_data.AddShallow("EMF", std::vector<std::string>{"B_CT.emf"}); // TODO this gets weird if we partition
                auto t_emf_local = tl.AddTask(t_flux_bounds, B_CT::CalculateEMF, md_sub_step_init.get());
                t_emf = KHARMADriver::AddBoundarySync(t_emf_local, tl, md_emf_only);
            }
            auto t_load_send_flux = tl.AddTask(t_emf, parthenon::LoadAndSendFluxCorrections, md_sub_step_init);
            auto t_recv_flux = tl.AddTask(t_load_send_flux, parthenon::ReceiveFluxCorrections, md_sub_step_init);
            t_flux_bounds = tl.AddTask(t_recv_flux, parthenon::SetFluxCorrections, md_sub_step_init);
        }

        // Apply the fluxes to calculate a change in cell-centered values "md_flux_src"
        auto t_flux_div = tl.AddTask(t_flux_bounds, Update::FluxDivergence<MeshData<Real>>, md_sub_step_init.get(), md_flux_src.get());

        // Add any source terms: geometric \Gamma * T, wind, damping, etc etc
        auto t_sources = tl.AddTask(t_flux_div, Packages::AddSource, md_sub_step_init.get(), md_flux_src.get());

        // UPDATE VARIABLES
        // TODO abstract this since the drivers share it
        // This block is designed to intelligently update a set of variables partially marked "Implicit"
        // and partially "Explicit," by first doing any explicit updates, then using them as elements
        // of the "guess" for the implicit solve

        // Update the explicitly-evolved variables using the source term
        // Add any proportion of the step start required by the integrator (e.g., RK2)
        auto t_avg_data_c = tl.AddTask(t_sources, Update::WeightedSumData<std::vector<MetadataFlag>, MeshData<Real>>,
                                    std::vector<MetadataFlag>({Metadata::GetUserFlag("Explicit"), Metadata::Independent, Metadata::Cell}),
                                    md_sub_step_init.get(), md_full_step_init.get(),
                                    integrator->gam0[stage-1], integrator->gam1[stage-1],
                                    md_solver.get());
        auto t_avg_data = t_avg_data_c;
        if (use_b_ct) {
            t_avg_data = tl.AddTask(t_avg_data_c, WeightedSumDataFace,
                                    std::vector<MetadataFlag>({Metadata::GetUserFlag("Explicit"), Metadata::Independent, Metadata::Face}),
                                    md_sub_step_init.get(), md_full_step_init.get(),
                                    integrator->gam0[stage-1], integrator->gam1[stage-1],
                                    md_solver.get());
        }
        // apply du/dt to the result
        auto t_update_c = tl.AddTask(t_sources, Update::WeightedSumData<std::vector<MetadataFlag>, MeshData<Real>>,
                                    std::vector<MetadataFlag>({Metadata::GetUserFlag("Explicit"), Metadata::Independent, Metadata::Cell}),
                                    md_solver.get(), md_flux_src.get(),
                                    1.0, integrator->beta[stage-1] * integrator->dt,
                                    md_solver.get());
        auto t_update = t_update_c;
        if (use_b_ct) {
            t_update = tl.AddTask(t_update_c, WeightedSumDataFace,
                                  std::vector<MetadataFlag>({Metadata::GetUserFlag("Explicit"), Metadata::Independent, Metadata::Face}),
                                  md_solver.get(), md_flux_src.get(),
                                  1.0, integrator->beta[stage-1] * integrator->dt,
                                  md_solver.get());
        }

        // If evolving GRMHD explicitly, UtoP needs a guess in order to converge, so we copy in md_sub_step_init
        auto t_copy_prims = t_update;
        if (!pkgs.at("GRMHD")->Param<bool>("implicit")) {
            t_copy_prims = tl.AddTask(t_none, Copy<MeshData<Real>>, std::vector<MetadataFlag>({Metadata::GetUserFlag("HD"), Metadata::GetUserFlag("Primitive")}),
                                      md_sub_step_init.get(), md_solver.get());
        }

        // Make sure the primitive values of *explicitly-evolved* variables are updated.
        // Packages with implicitly-evolved vars should only register BoundaryUtoP or BoundaryPtoU
        auto t_explicit_UtoP = tl.AddTask(t_copy_prims | t_update, Packages::MeshUtoP, md_solver.get(), IndexDomain::entire, false);

        // Done with explicit update
        auto t_explicit = t_explicit_UtoP;

        auto t_implicit = t_explicit;
        if (use_implicit) {
            // Extra containers for implicit solve
            std::shared_ptr<MeshData<Real>> &md_linesearch = (use_linesearch) ? pmesh->mesh_data.GetOrAdd("linesearch", i) : md_solver;

            // Copy the current state of any implicitly-evolved vars (at least the prims) in as a guess.
            // This sets md_solver = md_sub_step_init
            auto t_copy_guess = tl.AddTask(t_sources, Copy<MeshData<Real>>, std::vector<MetadataFlag>({Metadata::GetUserFlag("Implicit")}),
                                        md_sub_step_init.get(), md_solver.get());

            auto t_guess_ready = t_explicit | t_copy_guess;

            // The `solver` MeshData object now has the implicit primitives corresponding to initial/half step and
            // explicit variables have been updated to match the current step.
            // Copy the primitives to the `linesearch` MeshData object if linesearch was enabled.
            auto t_copy_linesearch = t_guess_ready;
            if (use_linesearch) {
                t_copy_linesearch = tl.AddTask(t_guess_ready, Copy<MeshData<Real>>, std::vector<MetadataFlag>({Metadata::GetUserFlag("Primitive")}),
                                                md_solver.get(), md_linesearch.get());
            }


            // Time-step implicit variables by root-finding the residual.
            // This calculates the primitive values after the substep for all "isImplicit" variables --
            // no need for separately adding the flux divergence or calling UtoP
            auto t_implicit_step = tl.AddTask(t_copy_linesearch, Implicit::Step, md_full_step_init.get(), md_sub_step_init.get(), 
                                         md_flux_src.get(), md_linesearch.get(), md_solver.get(), integrator->beta[stage-1] * integrator->dt);

            // Copy the entire solver state (everything defined on the grid, incl. our new Face variables) into the final state md_sub_step_final
            // If we're entirely explicit, we just declare these equal
            auto t_implicit_c = tl.AddTask(t_implicit_step, Copy<MeshData<Real>>, std::vector<MetadataFlag>({Metadata::Cell}),
                                    md_solver.get(), md_sub_step_final.get());
            t_implicit = tl.AddTask(t_implicit_step, WeightedSumDataFace, std::vector<MetadataFlag>({Metadata::Face}),
                                    md_solver.get(), md_solver.get(), 1.0, 0.0, md_sub_step_final.get());
        }

        // Apply all floors & limits (GRMHD,EMHD,etc), but do *not* immediately correct UtoP failures with FixUtoP --
        // rather, we will synchronize (including pflags!) first.
        // With an extra ghost zone, this *should* still allow binary-similar evolution between numbers of mesh blocks,
        // but hasn't been tested to do so yet.
        auto t_floors = tl.AddTask(t_implicit, Packages::MeshApplyFloors, md_sub_step_final.get(), IndexDomain::interior);

        KHARMADriver::AddBoundarySync(t_floors, tl, md_sync);
    }

    // Async Region: Any post-sync tasks.  Fixups, timestep & AMR tagging.
    //TaskRegion &async_region2 = tc.AddRegion(blocks.size());
    TaskRegion &async_region2 = tc.AddRegion(num_partitions);
    for (int i = 0; i < num_partitions; i++) {
        auto &tl  = async_region2[i];
        auto &md_sub_step_init  = pmesh->mesh_data.GetOrAdd(integrator->stage_name[stage-1], i);
        auto &md_sub_step_final = pmesh->mesh_data.GetOrAdd(integrator->stage_name[stage], i);
        auto &md_sync = pmesh->mesh_data.AddShallow("sync"+integrator->stage_name[stage]+std::to_string(i), md_sub_step_final, sync_vars);

        // If we're evolving the GRMHD variables explicitly, we need to fix UtoP variable inversion failures.
        // If implicitly, we run a (very similar) fix for solver failures.
        // Syncing bounds before calling this, and then running it over the whole domain, will make
        // behavior for different mesh breakdowns much more similar (identical?), since bad zones in
        // relevant ghost zone ranks will get to use all the same neighbors as if they were in the bulk
        // TODO fixups as a callback?
        auto t_fix_utop = t_none;
        if (!pkgs.at("GRMHD")->Param<bool>("implicit")) {
            t_fix_utop = tl.AddTask(t_none, Inverter::MeshFixUtoP, md_sub_step_final.get());
        }
        auto t_fix_solve = t_fix_utop;
        if (use_implicit) {
            t_fix_solve = tl.AddTask(t_fix_utop, Implicit::MeshFixSolve, md_sub_step_final.get());
        }

        // Re-apply boundary conditions to reflect fixes
        auto t_set_bc = tl.AddTask(t_fix_solve, parthenon::ApplyBoundaryConditionsOnCoarseOrFineMD, md_sync, false);

        // Any package- (likely, problem-) specific source terms which must be applied to primitive variables
        // Apply these only after the final step so they're operator-split
        auto t_prim_source = t_set_bc;
        if (stage == integrator->nstages) {
            t_prim_source = tl.AddTask(t_set_bc, Packages::MeshApplyPrimSource, md_sub_step_final.get());
        }
        // Electron heating goes where it does in the KHARMA Driver, for the same reasons
        auto t_heat_electrons = t_prim_source;
        if (use_electrons) {
            t_heat_electrons = tl.AddTask(t_prim_source, Electrons::MeshApplyElectronHeating,
                                          md_sub_step_init.get(), md_sub_step_final.get());
        }

        // Make sure *all* conserved vars are synchronized at step end
        auto t_ptou = tl.AddTask(t_heat_electrons, Flux::MeshPtoU, md_sub_step_final.get(), IndexDomain::entire, false);

        auto t_step_done = t_ptou;

        // Estimate next time step based on ctop
        if (stage == integrator->nstages) {
            auto t_new_dt =
                tl.AddTask(t_step_done, Update::EstimateTimestep<MeshData<Real>>, md_sub_step_final.get());

            // Update refinement
            if (pmesh->adaptive) {
                auto tag_refine = tl.AddTask(
                    t_step_done, parthenon::Refinement::Tag<MeshData<Real>>, md_sub_step_final.get());
            }
        }
    }

    // B Field cleanup: this is a separate solve so it's split out
    // It's also really slow when enabled so we don't care too much about limiting regions, etc.
    if (use_b_cleanup && (stage == integrator->nstages) && B_Cleanup::CleanupThisStep(pmesh, tm.ncycle)) {
        TaskRegion &cleanup_region = tc.AddRegion(num_partitions);
        for (int i = 0; i < num_partitions; i++) {
            auto &tl = cleanup_region[i];
            auto &md_sub_step_final = pmesh->mesh_data.GetOrAdd(integrator->stage_name[stage], i);
            tl.AddTask(t_none, B_Cleanup::CleanupDivergence, md_sub_step_final);
        }
    }

    // Second boundary sync:
    // ensure that primitive variables in ghost zones are *exactly*
    // identical to their physical counterparts, now that they have been
    // modified on each rank.
    const auto &two_sync = pkgs.at("Driver")->Param<bool>("two_sync");
    if (two_sync) {
        TaskRegion &async_region3 = tc.AddRegion(num_partitions);
        for (int i = 0; i < num_partitions; i++) {
            auto &md_sub_step_final = pmesh->mesh_data.GetOrAdd(integrator->stage_name[stage], i);
            auto &md_sync = pmesh->mesh_data.AddShallow("sync"+integrator->stage_name[stage]+std::to_string(i), md_sub_step_final, sync_vars);
            KHARMADriver::AddFullSyncRegion(tc, md_sync);
        }
    }

    return tc;
}

