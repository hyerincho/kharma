/* 
 *  File: b_flux_ct.cpp
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
#include "b_ct.hpp"

#include "decs.hpp"
#include "domain.hpp"
#include "grmhd.hpp"
#include "kharma.hpp"
// TODO eliminate sync
#include "kharma_driver.hpp"

#include <parthenon/parthenon.hpp>
#include <prolong_restrict/pr_ops.hpp>

using namespace parthenon;
using parthenon::refinement_ops::ProlongateSharedMinMod;
using parthenon::refinement_ops::RestrictAverage;
using parthenon::refinement_ops::ProlongateInternalAverage;

std::shared_ptr<KHARMAPackage> B_CT::Initialize(ParameterInput *pin, std::shared_ptr<Packages_t>& packages)
{
    auto pkg = std::make_shared<KHARMAPackage>("B_CT");
    Params &params = pkg->AllParams();

    // Diagnostic & inadvisable flags

    // KHARMA requires some kind of field transport if there is a magnetic field allocated.
    // Use this flag if you actually want to disable all magnetic field flux corrections,
    // and allow a field divergence to grow unchecked, usually for debugging or comparison reasons
    bool disable_ct = pin->GetOrAddBoolean("b_field", "disable_ct", false);
    params.Add("disable_ct", disable_ct);

    // Default to stopping execution when divB is large, which generally indicates something
    // has gone wrong.  As always, can be disabled by the brave.
    bool kill_on_large_divb = pin->GetOrAddBoolean("b_field", "kill_on_large_divb", true);
    params.Add("kill_on_large_divb", kill_on_large_divb);
    Real kill_on_divb_over = pin->GetOrAddReal("b_field", "kill_on_divb_over", 1.e-3);
    params.Add("kill_on_divb_over", kill_on_divb_over);

    // Currently bs99, gs05_c, gs05_0
    // TODO gs05_alpha, LDZ04 UCT1, LDZ07 UCT2
    std::string ct_scheme = pin->GetOrAddString("b_field", "ct_scheme", "bs99");
    params.Add("ct_scheme", ct_scheme);
    // Use the default Parthenon prolongation operator, rather than the divergence-preserving one
    // This relies entirely on the EMF communication for preserving the divergence
    bool lazy_prolongation = pin->GetOrAddBoolean("b_field", "lazy_prolongation", false);
    // Need to preserve divergence if you refine/derefine during sim i.e. AMR
    if (lazy_prolongation && pin->GetString("parthenon/mesh", "refinement") == "adaptive")
        throw std::runtime_error("Cannot use non-preserving prolongation in AMR!");
    
    // FIELDS

    // Flags for B fields on faces.
    // We don't mark these as "Conserved" else they'd be bundled
    // with all the cell vars in a bunch of places we don't want
    // Also note we *always* sync B field conserved var
    std::vector<MetadataFlag> flags_cons_f = {Metadata::Real, Metadata::Face, Metadata::Independent, Metadata::Restart,
                                              Metadata::GetUserFlag("Explicit"), Metadata::FillGhost};
    auto m = Metadata(flags_cons_f);
    if (!lazy_prolongation)
        m.RegisterRefinementOps<ProlongateSharedMinMod, RestrictAverage, ProlongateInternalOlivares>();
    else
        m.RegisterRefinementOps<ProlongateSharedMinMod, RestrictAverage, ProlongateInternalAverage>();
    pkg->AddField("cons.fB", m);

    // Cell-centered versions.  Needed for BS, not for other schemes.
    // Probably will want to keep primitives for e.g. correct PtoU of MHD vars, but cons maybe can go
    std::vector<MetadataFlag> flags_prim = {Metadata::Real, Metadata::Cell, Metadata::Derived, Metadata::GetUserFlag("Primitive"),
                                            Metadata::GetUserFlag("MHD"), Metadata::GetUserFlag("Explicit"), Metadata::Vector, Metadata::Restart};
    std::vector<MetadataFlag> flags_cons = {Metadata::Real, Metadata::Cell, Metadata::Derived, Metadata::Conserved, Metadata::WithFluxes,
                                            Metadata::GetUserFlag("MHD"), Metadata::GetUserFlag("Explicit"), Metadata::Vector};
    std::vector<int> s_vector({NVEC});
    m = Metadata(flags_prim, s_vector);
    pkg->AddField("prims.B", m);
    m = Metadata(flags_cons, s_vector);
    pkg->AddField("cons.B", m);
    if (packages->Get("Globals")->Param<std::string>("problem") == "resize_restart_kharma") {
        m = Metadata({Metadata::Real, Metadata::Face, Metadata::Derived, Metadata::FillGhost});
        pkg->AddField("B_Save", m);
    }

    // EMF on edges.
    std::vector<MetadataFlag> flags_emf = {Metadata::Real, Metadata::Edge, Metadata::Derived, Metadata::OneCopy, Metadata::FillGhost};
    m = Metadata(flags_emf);
    pkg->AddField("B_CT.emf", m);

    if (ct_scheme != "bs99") {
        std::vector<MetadataFlag> flags_emf_c = {Metadata::Real, Metadata::Cell, Metadata::Derived, Metadata::OneCopy};
        m = Metadata(flags_emf_c, s_vector);
        pkg->AddField("B_CT.cemf", m);
    }

    // Hyerin (04/04/24) averaged B fields needed for ismr
    m = Metadata({Metadata::Real, Metadata::Face, Metadata::Derived, Metadata::FillGhost});
    pkg->AddField("fB_avg", m);

    // CALLBACKS

    // We implement a source term replacement, rather than addition,
    // but same difference really
    pkg->AddSource = B_CT::AddSource;

    // Also ensure that prims get filled, both during step and on boundaries
    //pkg->MeshUtoP = B_CT::MeshUtoP;
    pkg->BlockUtoP = B_CT::BlockUtoP;
    pkg->BoundaryUtoP = B_CT::BlockUtoP;

    // Register the other callbacks
    pkg->PostStepDiagnosticsMesh = B_CT::PostStepDiagnostics;

    // The definition of MaxDivB we care about actually changes per-transport,
    // so calculating it is handled by the transport package
    // We'd only ever need to declare or calculate divB for output (getting the max is independent)
    if (KHARMA::FieldIsOutput(pin, "divB")) {
        pkg->BlockUserWorkBeforeOutput = B_CT::FillOutput;
        m = Metadata({Metadata::Real, Metadata::Cell, Metadata::Derived, Metadata::OneCopy});
        pkg->AddField("divB", m);
    }

    // List (vector) of HistoryOutputVars that will all be enrolled as output variables
    // LATER
    parthenon::HstVar_list hst_vars = {};
    hst_vars.emplace_back(parthenon::HistoryOutputVar(UserHistoryOperation::max, B_CT::MaxDivB, "MaxDivB"));
    // Event horizon magnetization.  Might be the same or different for different representations?
    if (pin->GetBoolean("coordinates", "spherical")) {
        // hst_vars.emplace_back(parthenon::HistoryOutputVar(UserHistoryOperation::sum, ReducePhi0, "Phi_0"));
        // hst_vars.emplace_back(parthenon::HistoryOutputVar(UserHistoryOperation::sum, ReducePhi5, "Phi_EH"));
    }
    // add callbacks for HST output to the Params struct, identified by the `hist_param_key`
    pkg->AddParam<>(parthenon::hist_param_key, hst_vars);

    return pkg;
}

TaskStatus B_CT::MeshUtoP(MeshData<Real> *md, IndexDomain domain, bool coarse)
{
    // TODO later
    for (int i=0; i < md->NumBlocks(); i++)
        B_CT::BlockUtoP(md->GetBlockData(i).get(), domain, coarse);
    return TaskStatus::complete;
}

TaskStatus B_CT::BlockUtoP(MeshBlockData<Real> *rc, IndexDomain domain, bool coarse)
{
    auto pmb = rc->GetBlockPointer();
    const int ndim = pmb->pmy_mesh->ndim;
    auto B_Uf = rc->PackVariables(std::vector<std::string>{"cons.fB"});
    auto B_U = rc->PackVariables(std::vector<std::string>{"cons.B"});
    auto B_P = rc->PackVariables(std::vector<std::string>{"prims.B"});
    const auto& G = pmb->coords;
    // Return if we're not syncing U & P at all (e.g. edges)
    if (B_Uf.GetDim(4) == 0) return TaskStatus::complete;

    const IndexRange3 bc = KDomain::GetRange(rc, domain, coarse);

    // Average the primitive vals to zone centers
    pmb->par_for("UtoP_B_center", bc.ks, bc.ke, bc.js, bc.je, bc.is, bc.ie,
        KOKKOS_LAMBDA (const int &k, const int &j, const int &i) {
            B_P(V1, k, j, i) = (B_Uf(F1, 0, k, j, i) / G.gdet(Loci::face1, j, i)
                              + B_Uf(F1, 0, k, j, i + 1) / G.gdet(Loci::face1, j, i + 1)) / 2;
            B_P(V2, k, j, i) = (ndim > 1) ? (B_Uf(F2, 0, k, j, i) / G.gdet(Loci::face2, j, i)
                                           + B_Uf(F2, 0, k, j + 1, i) / G.gdet(Loci::face2, j + 1, i)) / 2
                                           : B_Uf(F2, 0, k, j, i) / G.gdet(Loci::face2, j, i);
            B_P(V3, k, j, i) = (ndim > 2) ? (B_Uf(F3, 0, k, j, i) / G.gdet(Loci::face3, j, i)
                                           + B_Uf(F3, 0, k + 1, j, i) / G.gdet(Loci::face3, j, i)) / 2
                                          : B_Uf(F3, 0, k, j, i) / G.gdet(Loci::face3, j, i);
        }
    );
    // Recover conserved B at centers
    pmb->par_for("UtoP_B_centerPtoU", 0, NVEC-1, bc.ks, bc.ke, bc.js, bc.je, bc.is, bc.ie,
        KOKKOS_LAMBDA (const int &v, const int &k, const int &j, const int &i) {
            B_U(v, k, j, i) = B_P(v, k, j, i) * G.gdet(Loci::center, j, i);
        }
    );

    return TaskStatus::complete;
}

TaskStatus B_CT::CalculateEMF(MeshData<Real> *md)
{
    auto pmesh = md->GetMeshPointer();
    const int ndim = pmesh->ndim;

    // EMF temporary
    auto& emf_pack = md->PackVariables(std::vector<std::string>{"B_CT.emf"});

    // Figure out indices
    const IndexRange3 b = KDomain::GetRange(md, IndexDomain::interior, 0, 0);
    const IndexRange3 b1 = KDomain::GetRange(md, IndexDomain::interior, 0, 1);
    const IndexRange block = IndexRange{0, emf_pack.GetDim(5)-1};

    auto pmb0 = md->GetBlockData(0)->GetBlockPointer();

    // Calculate circulation by averaging fluxes
    // This is the base of most other schemes, which make corrections
    // It is the entirety of B&S '99
    auto& B_U = md->PackVariablesAndFluxes(std::vector<std::string>{"cons.B"});
    pmb0->par_for("B_CT_emf_BS", block.s, block.e, b1.ks, b1.ke, b1.js, b1.je, b1.is, b1.ie,
        KOKKOS_LAMBDA (const int &bl, const int &k, const int &j, const int &i) {
            // The basic EMF per length along edges is the B field flux
            // We use this form rather than multiply by edge length here,
            // since the default restriction op averages values
            const auto& G = B_U.GetCoords(bl);
            if (ndim > 2) {
                emf_pack(bl, E1, 0, k, j, i) =
                    0.25*(B_U(bl).flux(X2DIR, V3, k - 1, j, i) + B_U(bl).flux(X2DIR, V3, k, j, i)
                        - B_U(bl).flux(X3DIR, V2, k, j - 1, i) - B_U(bl).flux(X3DIR, V2, k, j, i));
                emf_pack(bl, E2, 0, k, j, i) =
                    0.25*(B_U(bl).flux(X3DIR, V1, k, j, i - 1) + B_U(bl).flux(X3DIR, V1, k, j, i)
                        - B_U(bl).flux(X1DIR, V3, k - 1, j, i) - B_U(bl).flux(X1DIR, V3, k, j, i));
                emf_pack(bl, E3, 0, k, j, i) =
                    0.25*(B_U(bl).flux(X1DIR, V2, k, j - 1, i) + B_U(bl).flux(X1DIR, V2, k, j, i)
                        - B_U(bl).flux(X2DIR, V1, k, j, i - 1) - B_U(bl).flux(X2DIR, V1, k, j, i));
            } else if (ndim > 1) {
                emf_pack(bl, E1, 0, k, j, i) =  B_U(bl).flux(X2DIR, V3, k, j, i);
                emf_pack(bl, E2, 0, k, j, i) = -B_U(bl).flux(X1DIR, V3, k, j, i);
                emf_pack(bl, E3, 0, k, j, i) =
                    0.25*(B_U(bl).flux(X1DIR, V2, k, j - 1, i) + B_U(bl).flux(X1DIR, V2, k, j, i)
                        - B_U(bl).flux(X2DIR, V1, k, j, i - 1) - B_U(bl).flux(X2DIR, V1, k, j, i));
            } else {
                emf_pack(bl, E1, 0, k, j, i) = 0;
                emf_pack(bl, E2, 0, k, j, i) = -B_U(bl).flux(X1DIR, V3, k, j, i);
                emf_pack(bl, E3, 0, k, j, i) =  B_U(bl).flux(X1DIR, V2, k, j, i);
            }
        }
    );

    std::string scheme = pmesh->packages.Get("B_CT")->Param<std::string>("ct_scheme");
    if (scheme == "bs99") {
        // Nothing more to do
    } else if (scheme == "gs05_0" || scheme == "gs05_c") {
        // Additional terms for Stone & Gardiner '09
        // Average fluxes and derivatives
        auto& uvec = md->PackVariables(std::vector<std::string>{"prims.uvec"});
        auto& emfc = md->PackVariables(std::vector<std::string>{"B_CT.cemf"});
        auto& B_U = md->PackVariablesAndFluxes(std::vector<std::string>{"cons.B"});
        auto& B_P = md->PackVariables(std::vector<std::string>{"prims.B"});
        // emf in center == -v x B
        pmb0->par_for("B_CT_emfc", block.s, block.e, b.ks, b.ke, b.js, b.je, b.is, b.ie,
            KOKKOS_LAMBDA (const int &bl, const int &k, const int &j, const int &i) {
                VLOOP emfc(bl, v, k, j, i) = 0.;
                VLOOP3 emfc(bl, x, k, j, i) -= antisym(v, w, x) * uvec(bl, v, k, j, i) * B_U(bl, w, k, j, i);
            }
        );

        if (scheme == "gs05_0") {
            const int kd = ndim > 2 ? 1 : 0;
            const int jd = ndim > 1 ? 1 : 0;
            const int id = ndim > 0 ? 1 : 0;
            pmb0->par_for("B_CT_emf_GS05_0", block.s, block.e, b1.ks, b1.ke, b1.js, b1.je, b1.is, b1.ie,
                KOKKOS_LAMBDA (const int &bl, const int &k, const int &j, const int &i) {
                    const auto& G = B_U.GetCoords(bl);
                    // Just subtract centered emf from twice the face version
                    // More stable for planar flows even without anything fancy
                    emf_pack(bl, E1, 0, k, j, i) = 2 * emf_pack(bl, E1, 0, k, j, i)
                        - 0.25*(emfc(bl, V1, k, j, i)      + emfc(bl, V1, k, j - jd, i)
                              + emfc(bl, V1, k, j - jd, i) + emfc(bl, V1, k - kd, j - jd, i));
                    emf_pack(bl, E2, 0, k, j, i) = 2 * emf_pack(bl, E2, 0, k, j, i)
                        - 0.25*(emfc(bl, V2, k, j, i)      + emfc(bl, V2, k, j, i - id)
                              + emfc(bl, V2, k - kd, j, i) + emfc(bl, V2, k - kd, j, i - id));
                    emf_pack(bl, E3, 0, k, j, i) = 2 * emf_pack(bl, E3, 0, k, j, i)
                        - 0.25*(emfc(bl, V3, k, j, i)      + emfc(bl, V3, k, j, i - id)
                              + emfc(bl, V3, k, j - jd, i) + emfc(bl, V3, k, j - jd, i - id));
                }
            );
        } else if (scheme == "gs05_c") {
            // Get primitive velocity at face (on right side) (TODO do we need some average?)
            auto& uvecf = md->PackVariables(std::vector<std::string>{"Flux.vr"});

            pmb0->par_for("B_CT_emf_GS05_c", block.s, block.e, b1.ks, b1.ke, b1.js, b1.je, b1.is, b1.ie,
                KOKKOS_LAMBDA (const int &bl, const int &k, const int &j, const int &i) {
                    const auto& G = B_U.GetCoords(bl);
                    // "simple" flux + upwinding method, Stone & Gardiner '09 but also in Stone+08 etc.
                    // Upwinded differences take in order (1-indexed):
                    // 1. EMF component direction to calculate
                    // 2. Direction of derivative
                    // 3. Direction of upwinding
                    // ...then zone number...
                    // and finally, a boolean indicating a leftward (e.g., i-3/4) vs rightward (i-1/4) position
                    // TODO(BSP) This doesn't properly support 2D. Yell when it's chosen?
                    if (ndim > 2) {
                        emf_pack(bl, E1, 0, k, j, i) +=
                              0.25*(upwind_diff(B_U(bl), emfc(bl), uvecf(bl), 1, 3, 2, k, j, i, false)
                                  - upwind_diff(B_U(bl), emfc(bl), uvecf(bl), 1, 3, 2, k, j, i, true))
                            + 0.25*(upwind_diff(B_U(bl), emfc(bl), uvecf(bl), 1, 2, 3, k, j, i, false)
                                  - upwind_diff(B_U(bl), emfc(bl), uvecf(bl), 1, 2, 3, k, j, i, true));
                        emf_pack(bl, E2, 0, k, j, i) +=
                              0.25*(upwind_diff(B_U(bl), emfc(bl), uvecf(bl), 2, 1, 3, k, j, i, false)
                                  - upwind_diff(B_U(bl), emfc(bl), uvecf(bl), 2, 1, 3, k, j, i, true))
                            + 0.25*(upwind_diff(B_U(bl), emfc(bl), uvecf(bl), 2, 3, 1, k, j, i, false)
                                  - upwind_diff(B_U(bl), emfc(bl), uvecf(bl), 2, 3, 1, k, j, i, true));
                    }
                    emf_pack(bl, E3, 0, k, j, i) +=
                          0.25*(upwind_diff(B_U(bl), emfc(bl), uvecf(bl), 3, 2, 1, k, j, i, false)
                              - upwind_diff(B_U(bl), emfc(bl), uvecf(bl), 3, 2, 1, k, j, i, true))
                        + 0.25*(upwind_diff(B_U(bl), emfc(bl), uvecf(bl), 3, 1, 2, k, j, i, false)
                              - upwind_diff(B_U(bl), emfc(bl), uvecf(bl), 3, 1, 2, k, j, i, true));
                }
            );
        }
    } else {
        throw std::invalid_argument("Invalid CT scheme specified!  Must be one of bs99, gs05_0, gs05_c!");
    }

    return TaskStatus::complete;
}

TaskStatus B_CT::AddSource(MeshData<Real> *md, MeshData<Real> *mdudt)
{
    auto pmesh = md->GetMeshPointer();
    const int ndim = pmesh->ndim;

    // EMF temporary
    auto& emf_pack = md->PackVariables(std::vector<std::string>{"B_CT.emf"});

    // Figure out indices
    const IndexRange3 b = KDomain::GetRange(md, IndexDomain::interior, 0, 0);
    const IndexRange3 b1 = KDomain::GetRange(md, IndexDomain::interior, 0, 1);
    const IndexRange block = IndexRange{0, emf_pack.GetDim(5)-1};

    auto pmb0 = md->GetBlockData(0)->GetBlockPointer();

    // This is what we're replacing
    auto& dB_Uf_dt = mdudt->PackVariables(std::vector<std::string>{"cons.fB"});
    // Circulation -> change in flux at face
    pmb0->par_for("B_CT_Circ_1", block.s, block.e, b.ks, b.ke, b.js, b.je, b1.is, b1.ie,
        KOKKOS_LAMBDA (const int &bl, const int &k, const int &j, const int &i) {
            const auto& G = dB_Uf_dt.GetCoords(bl);
            dB_Uf_dt(bl, F1, 0, k, j, i) = (G.Volume<E3>(k, j + 1, i) * emf_pack(bl, E3, 0, k, j + 1, i)
                                          - G.Volume<E3>(k, j, i)     * emf_pack(bl, E3, 0, k, j, i));
            if (ndim > 2)
                dB_Uf_dt(bl, F1, 0, k, j, i) += (-G.Volume<E2>(k + 1, j, i) * emf_pack(bl, E2, 0, k + 1, j, i)
                                                + G.Volume<E2>(k, j, i)     * emf_pack(bl, E2, 0, k, j, i));
            dB_Uf_dt(bl, F1, 0, k, j, i) /= G.Volume<F1>(k, j, i);
        }
    );
    pmb0->par_for("B_CT_Circ_2", block.s, block.e, b.ks, b.ke, b1.js, b1.je, b.is, b.ie,
        KOKKOS_LAMBDA (const int &bl, const int &k, const int &j, const int &i) {
            const auto& G = dB_Uf_dt.GetCoords(bl);
            dB_Uf_dt(bl, F2, 0, k, j, i) = (-G.Volume<E3>(k, j, i + 1) * emf_pack(bl, E3, 0, k, j, i + 1)
                                           + G.Volume<E3>(k, j, i)     * emf_pack(bl, E3, 0, k, j, i));
            if (ndim > 2)
                dB_Uf_dt(bl, F2, 0, k, j, i) +=  (G.Volume<E1>(k + 1, j, i) * emf_pack(bl, E1, 0, k + 1, j, i)
                                                - G.Volume<E1>(k, j, i)     * emf_pack(bl, E1, 0, k, j, i));
            dB_Uf_dt(bl, F2, 0, k, j, i) /= G.Volume<F2>(k, j, i);
        }
    );
    pmb0->par_for("B_CT_Circ_3", block.s, block.e, b1.ks, b1.ke, b.js, b.je, b.is, b.ie,
        KOKKOS_LAMBDA (const int &bl, const int &k, const int &j, const int &i) {
            const auto& G = dB_Uf_dt.GetCoords(bl);
            dB_Uf_dt(bl, F3, 0, k, j, i) = (G.Volume<E2>(k, j, i + 1) * emf_pack(bl, E2, 0, k, j, i + 1)
                                            - G.Volume<E2>(k, j, i)     * emf_pack(bl, E2, 0, k, j, i)
                                            - G.Volume<E1>(k, j + 1, i) * emf_pack(bl, E1, 0, k, j + 1, i)
                                            + G.Volume<E1>(k, j, i)     * emf_pack(bl, E1, 0, k, j, i)) / G.Volume<F3>(k, j, i);
        }
    );

    return TaskStatus::complete;
}

double B_CT::MaxDivB(MeshData<Real> *md)
{
    auto pmesh = md->GetMeshPointer();
    const int ndim = pmesh->ndim;

    auto B_U = md->PackVariables(std::vector<std::string>{"cons.fB"});

    // Figure out indices
    const IndexRange ib = md->GetBoundsI(IndexDomain::interior);
    const IndexRange jb = md->GetBoundsJ(IndexDomain::interior);
    const IndexRange kb = md->GetBoundsK(IndexDomain::interior);
    const IndexRange block = IndexRange{0, B_U.GetDim(5)-1};

    auto pmb0 = md->GetBlockData(0)->GetBlockPointer();

    double max_divb;
    Kokkos::Max<double> max_reducer(max_divb);
    pmb0->par_reduce("divB_max", block.s, block.e, kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
        KOKKOS_LAMBDA (const int &b, const int &k, const int &j, const int &i, double &local_result) {
            const auto& G = B_U.GetCoords(b);
            double local_divb = face_div(G, B_U(b), ndim, k, j, i);
            if (local_divb > local_result) local_result = local_divb;
        }
    , max_reducer);

    return max_divb;
}
double B_CT::BlockMaxDivB(MeshBlockData<Real> *rc)
{
    const int ndim = KDomain::GetNDim(rc);

    auto B_U = rc->PackVariables(std::vector<std::string>{"cons.fB"});

    // Figure out indices
    const IndexRange3 b = KDomain::GetRange(rc, IndexDomain::interior);

    auto pmb = rc->GetBlockPointer();

    double max_divb;
    Kokkos::Max<double> max_reducer(max_divb);
    pmb->par_reduce("divB_max", b.ks, b.ke, b.js, b.je, b.is, b.ie,
        KOKKOS_LAMBDA (const int &k, const int &j, const int &i, double &local_result) {
            const auto& G = B_U.GetCoords();
            double local_divb = face_div(G, B_U, ndim, k, j, i);
            if (local_divb > local_result) local_result = local_divb;
        }
    , max_reducer);

    return max_divb;
}

double B_CT::GlobalMaxDivB(MeshData<Real> *md, bool all_reduce)
{
    if (all_reduce) {
        Reductions::StartToAll<Real>(md, 2, MaxDivB(md), MPI_MAX);
        return Reductions::CheckOnAll<Real>(md, 2);
    } else {
        Reductions::Start<Real>(md, 2, MaxDivB(md), MPI_MAX);
        return Reductions::Check<Real>(md, 2);
    }
}

TaskStatus B_CT::PrintGlobalMaxDivB(MeshData<Real> *md, bool kill_on_large_divb)
{
    auto pmb0 = md->GetBlockData(0)->GetBlockPointer();

    // Since this is in the history file now, I don't bother printing it
    // unless we're being verbose. It's not costly to calculate though
    const bool print = pmb0->packages.Get("Globals")->Param<int>("verbose") >= 1;
    if (print || kill_on_large_divb) {
        // Calculate the maximum from/on all nodes
        const double divb_max = B_CT::GlobalMaxDivB(md);
        // Print on rank zero
        if (MPIRank0() && print) {
            printf("Max DivB: %g\n", divb_max); // someday I'll learn stream options
        }
        if (kill_on_large_divb) {
            if (divb_max > pmb0->packages.Get("B_CT")->Param<Real>("kill_on_divb_over"))
                throw std::runtime_error("DivB exceeds maximum! Quitting...");
        }
    }

    return TaskStatus::complete;
}

// TODO unify these by adding FillOutputMesh option

void B_CT::CalcDivB(MeshData<Real> *md, std::string divb_field_name)
{
    auto pmesh = md->GetMeshPointer();
    const int ndim = pmesh->ndim;

    // Packing out here avoids frequent per-mesh packs.  Do we need to?
    auto B_U = md->PackVariables(std::vector<std::string>{"cons.fB"});
    auto divB = md->PackVariables(std::vector<std::string>{divb_field_name});

    const IndexRange ib = md->GetBoundsI(IndexDomain::interior);
    const IndexRange jb = md->GetBoundsJ(IndexDomain::interior);
    const IndexRange kb = md->GetBoundsK(IndexDomain::interior);
    const IndexRange block = IndexRange{0, B_U.GetDim(5)-1};

    auto pmb0 = md->GetBlockData(0)->GetBlockPointer();

    // See MaxDivB for details
    pmb0->par_for("calc_divB", block.s, block.e, kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
        KOKKOS_LAMBDA (const int &b, const int &k, const int &j, const int &i) {
            const auto& G = B_U.GetCoords(b);
            divB(b, 0, k, j, i) = face_div(G, B_U(b), ndim, k, j, i);
        }
    );
}

void B_CT::FillOutput(MeshBlock *pmb, ParameterInput *pin)
{
    auto rc = pmb->meshblock_data.Get();
    const int ndim = pmb->pmy_mesh->ndim;
    if (ndim < 2) return;

    auto B_U = rc->PackVariables(std::vector<std::string>{"cons.fB"});
    auto divB = rc->PackVariables(std::vector<std::string>{"divB"});

    const IndexRange ib = rc->GetBoundsI(IndexDomain::interior);
    const IndexRange jb = rc->GetBoundsJ(IndexDomain::interior);
    const IndexRange kb = rc->GetBoundsK(IndexDomain::interior);
    const IndexRange block = IndexRange{0, B_U.GetDim(5)-1};

    pmb->par_for("divB_output", kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
        KOKKOS_LAMBDA (const int &k, const int &j, const int &i) {
            const auto& G = B_U.GetCoords();
            divB(0, k, j, i) = face_div(G, B_U, ndim, k, j, i);
        }
    );
}

TaskStatus B_CT::DerefinePoles(MeshData<Real> *md, uint nlevels)
{
    // HYERIN (01/17/24) this routine is not general yet and only applies to polar boundaries for now.
    auto B_U = md->PackVariables(std::vector<std::string>{"cons.fB"});
    auto B_avg = md->PackVariables(std::vector<std::string>{"fB_avg"});
    auto rho_U = md->PackVariables(std::vector<std::string>{"cons.rho"});
    auto u_U = md->PackVariables(std::vector<std::string>{"cons.u"});
    auto uvec_U = md->PackVariables(std::vector<std::string>{"cons.uvec"});
    const IndexRange block = IndexRange{0, B_U.GetDim(5)-1};
    auto pmb0 = md->GetBlockData(0)->GetBlockPointer();
    
    // Figure out indices
    IndexRange3 bCC, bF1, bF2, bF3;
    int j_f, offset, point_out, jps, jp_now;
    for (int i = 0; i < BOUNDARY_NFACES; i++) {
        BoundaryFace bface = (BoundaryFace)i;
        auto bname = KBoundaries::BoundaryName(bface);
        auto bdir = KBoundaries::BoundaryDirection(bface);
        auto domain = KBoundaries::BoundaryDomain(bface);
        auto binner = KBoundaries::BoundaryIsInner(bface);
        if (bdir == X2DIR) {
            // indices
            bCC = KDomain::GetRange(md, domain, CC);
            bF1 = KDomain::GetRange(md, domain, F1);
            bF2 = KDomain::GetRange(md, domain, F2, (binner) ? 0 : -1, (binner) ? 1 : 0, false); 
            bF3 = KDomain::GetRange(md, domain, F3);
            j_f = (binner) ? bF2.je : bF2.js; // last physical face
            jps = (binner) ? j_f + (nlevels - 1) : j_f - (nlevels - 1); // start of the lowest level of derefinement
            const IndexRange j_p = IndexRange{(binner) ? j_f : jps, (binner) ? jps : j_f};  // Range of x2 to be de-refined
            offset = (binner) ? 1 : -1; // offset to read the physical face values

            // F1 average
            pmb0->par_for("B_CT_derefine_poles_avg_F1", block.s, block.e, bCC.ks, bCC.ke, j_p.s, j_p.e, bF1.is, bF1.ie,
                KOKKOS_LAMBDA (const int &bl, const int &k, const int &j, const int &i) {
                    int j_c, coarse_cell_len, ktemp, k_fine, k_start;
                    Real avg;
                    const auto& G = B_U.GetCoords(bl);

                    coarse_cell_len = m::pow(2, ((binner) ? jps - j : j - jps) + 1);
                    avg = 0.;
                    j_c = j + ((binner) ? 0 : -1); // cell center
                    k_fine = k % coarse_cell_len; // this fine cell's k-index within the coarse cell
                    k_start = k - k_fine; // starting k-index of the coarse cell
                    
                    // average over fine cells within the coarse cell we're in
                    for (ktemp = 0; ktemp < coarse_cell_len; ++ktemp) 
                        avg += B_U(bl)(F1, 0, k_start + ktemp, j_c, i) * G.Volume<F1>(k_start + ktemp, j_c, i);
                    avg /= coarse_cell_len;

                    B_avg(bl)(F1, 0, k, j_c, i) = avg / G.Volume<F1>(k, j_c, i);
                }
            );
            // F2 average
            pmb0->par_for("B_CT_derefine_poles_avg_F2", block.s, block.e, bCC.ks, bCC.ke, j_p.s, j_p.e, bCC.is, bCC.ie,
                KOKKOS_LAMBDA (const int &bl, const int &k, const int &j, const int &i) {
                    int coarse_cell_len, ktemp, k_fine, k_start;
                    Real avg;
                    const auto& G = B_U.GetCoords(bl);

                    coarse_cell_len = m::pow(2, ((binner) ? jps - j : j - jps) + 1);
                    avg = 0.;
                    k_fine = k % coarse_cell_len; // fine cell's k index within the coarse cell
                    k_start = k - k_fine; // starting k-index of the coarse cell
                            
                    if (j == j_f) { // The fine cells have 0 fluxes through the physical-ghost boundaries.
                        B_avg(bl)(F2, 0, k, j, i) = 0.;
                    } else { // average the fine cells
                        avg = 0.;
                        for (ktemp = 0; ktemp < coarse_cell_len; ++ktemp)
                            avg += B_U(bl)(F2, 0, k_start + ktemp, j, i) * G.Volume<F2>(k_start + ktemp, j, i);
                        avg /= coarse_cell_len;
                        
                        B_avg(bl)(F2, 0, k, j, i) = avg / G.Volume<F2>(k, j, i);
                    }
                }
            );
            // F3 average
            pmb0->par_for("B_CT_derefine_poles_avg_F3", block.s, block.e, bF3.ks, bF3.ke, j_p.s, j_p.e, bCC.is, bCC.ie,
                KOKKOS_LAMBDA (const int &bl, const int &k, const int &j, const int &i) {
                    int j_c, coarse_cell_len, c_half, ktemp, k_fine, k_start, k_half, k_end, current_lv;
                    Real avg, B_start, B_center, B_end;
                    const auto& G = B_U.GetCoords(bl);
                    
                    current_lv = ((binner) ? jps - j : j - jps); // the current level of derefinement at given j
                    c_half = m::pow(2, current_lv); // half of the coarse cell's length
                    coarse_cell_len = 2 * c_half;
                    j_c = j + ((binner) ? 0 : -1); // cell center
                    k_fine = k % coarse_cell_len; // this fine cell's k-index within the coarse cell
                    k_start = k - k_fine; // starting k-index of the coarse cell
                    k_half = k_start + c_half;
                    k_end  = k_start + coarse_cell_len; // end k-index of the coarse cell
                    
                    if (k % coarse_cell_len == 0) { // at the faces of the coarse cells. Don't modify them.
                        B_avg(bl)(F3, 0, k, j_c, i) = B_U(bl)(F3, 0, k, j_c, i);
                    else {
                        // F3: The internal faces will take care of the divB=0. The two faces of the coarse cell will remain unchanged.
                        //// First calculate the very central internal face. In other words, deal with the highest level internal face first.
                        //// Sum of F2 fluxes in the left and right half of the coarse cell each. 
                        Real c_left_v = 0., c_right_v = 0.;
                        for (int ktemp = 0; ktemp < c_half; ++ktemp) {
                            c_left_v += B_U(bl)(F2, 0, k_half -1 - ktemp, j + offset, i) * G.Volume<F2>(k_half - 1 - ktemp, j + offset, i);
                            c_right_v += B_U(bl)(F2, 0, k_half   + ktemp, j + offset, i) * G.Volume<F2>(k_half     + ktemp, j + offset, i);
                        }
                        B_start = B_U(bl)(F3, 0, k_start, j_c, i) * G.Volume<F3>(k_start, j_c, i);
                        B_end   = B_U(bl)(F3, 0, k_end,   j_c, i) * G.Volume<F3>(k_end,   j_c, i);
                        B_center = (B_start + B_end + point_out * (c_right_v - c_left_v)) / 2.;
                        if (k == 8 && i ==8 && j == j_f) printf("HYERIN: F3 internal %.3g\n", B_center);
                        
                        if (k == k_half) { // if at the center, then store the calculated value.
                            B_avg(bl)(F3, 0, k, j_c, i) = B_center / G.Volume<F3>(k, j_c, i);
                        } else if (k < k_half) { // interpolate between B_start and B_center
                            B_avg(bl)(F3, 0, k, j_c, i) = ((c_half - k_fine) * B_start + k_fine * B_center) / (c_half * G.Volume<F3>(k, j_c, i));
                        } else if (k > k_half) { // interpolate between B_end and B_center
                            B_avg(bl)(F3, 0, k, j_c, i) = ((k_fine - c_half) * B_end + (coarse_cell_len - k_fine) * B_center) / (c_half * G.Volume<F3>(k, j_c, i));
                        }
                    }
                }
            );
            
            // F1 write
            pmb0->par_for("B_CT_derefine_poles_F1", block.s, block.e, bCC.ks, bCC.ke, j_p.s, j_p.e, bF1.is, bF1.ie,
                KOKKOS_LAMBDA (const int &bl, const int &k, const int &j, const int &i) {
                    int j_c;
                    j_c = j + ((binner) ? 0 : -1); // cell center
                    B_U(bl)(F1, 0, k, j_c, i) = B_avg(bl)(F1, 0, k, j_c, i);
                }
            );
            // F2 write
            pmb0->par_for("B_CT_derefine_poles_F2", block.s, block.e, bCC.ks, bCC.ke, j_p.s, j_p.e, bCC.is, bCC.ie,
                KOKKOS_LAMBDA (const int &bl, const int &k, const int &j, const int &i) {
                    B_U(bl)(F2, 0, k, j, i) = B_avg(bl)(F2, 0, k, j, i);
                }
            );
            // F3 write
            pmb0->par_for("B_CT_derefine_poles_F3", block.s, block.e, bF3.ks, bF3.ke, j_p.s, j_p.e, bCC.is, bCC.ie,
                KOKKOS_LAMBDA (const int &bl, const int &k, const int &j, const int &i) {
                    int j_c;
                    j_c = j + ((binner) ? 0 : -1); // cell center
                    B_U(bl)(F3, 0, k, j_c, i) = B_avg(bl)(F3, 0, k, j_c, i);
                }
            );
            //bF2 = KDomain::GetRange(md, domain, F2, (binner) ? 0 : -1, (binner) ? 1 : 0, false); 
            //bF1 = KDomain::GetRange(md, domain, F1);
            //j_f = (binner) ? bF2.je : bF2.js; // last physical face
            //jps = (binner) ? j_f + (nlevels - 1) : j_f - (nlevels - 1); // start of the lowest level of derefinement
            ////auto j_c = (binner) ? j_f : j_f - 1; // last physical cell
            //offset = (binner) ? 1 : -1; // offset to read the physical face values
            //point_out = offset; // if F2 B field at j_f + offset face is positive when pointing out of the cell, +1.
            //// layers near the poles to be derefined
            ////int jpe = j_f;
            //const IndexRange j_p = IndexRange{(binner) ? j_f : jps, (binner) ? jps : j_f}; 
            //for (int jtemp = 0; jtemp < nlevels; jtemp++) {
            //    jp_now = jps - point_out * jtemp; // The j index of the polar region that we are working on now.
            //    pmb0->par_for("B_CT_derefine_poles", block.s, block.e, bF1.ks, bF1.ke, jp_now, jp_now, bF1.is, bF1.ie,
            //        KOKKOS_LAMBDA (const int &bl, const int &k, const int &j, const int &i) {
            //            int j_c, coarse_cell_len, c_half, ktemp;
            //            Real avg;
            //            const auto& G = B_U.GetCoords(bl);
            //            
            //            coarse_cell_len = 2 * (((binner) ? jps - j : j - jps) + 1);
            //            c_half = coarse_cell_len / 2; // half of the coarse cell's length
            //            if (k % coarse_cell_len == 0) {
            //                // F1: just average over the fine cells
            //                avg = 0.;
            //                j_c = j + (binner ? 0 : -1); // cell center
            //                for (ktemp = 0; ktemp < coarse_cell_len; ++ktemp) 
            //                    avg += B_U(bl)(F1, 0, k + ktemp, j_c, i) * G.Volume<F1>(k + ktemp, j_c, i);
            //                avg /= coarse_cell_len;
            //                for (ktemp = 0; ktemp < coarse_cell_len; ++ktemp)
            //                    B_U(bl)(F1, 0, k + ktemp, j_c, i) = avg / G.Volume<F1>(k + ktemp, j_c, i);
            //                
            //                // F2: 
            //                if (j == j_f) { //The two fine cells have 0 fluxes through the physical-ghost boundaries.
            //                    for (ktemp = 0; ktemp < coarse_cell_len; ++ktemp) B_U(bl)(F2, 0, k + ktemp, j, i) = 0.;
            //                } else { // average the cells
            //                    avg = 0.;
            //                    for (ktemp = 0; ktemp < coarse_cell_len; ++ktemp)
            //                        avg += B_U(bl)(F2, 0, k + ktemp, j, i) * G.Volume<F2>(k + ktemp, j, i);
            //                    avg /= coarse_cell_len;
            //                    for (ktemp = 0; ktemp < coarse_cell_len; ++ktemp)
            //                        B_U(bl)(F2, 0, k + ktemp, j, i) = avg / G.Volume<F2>(k + ktemp, j, i);
            //                }
            //                
            //                // F3: The internal face between the two fine cells will take care of the divB=0. The two other faces of the coarse cell will remain unchanged.
            //                // First modify the very central internal face. In other words, deal with the highest level internal face first.
            //                // Sum of F2 fluxes in the left and right half of the coarse cell each. 
            //                Real c_left_v = 0., c_right_v = 0.;
            //                for (int ktemp = 0; ktemp < c_half; ++ktemp) {
            //                    c_left_v += B_U(bl)(F2, 0, k + c_half - ktemp - 1, j + offset, i) * G.Volume<F2>(k + c_half - ktemp - 1, j + offset, i);
            //                    c_right_v += B_U(bl)(F2, 0, k + c_half + ktemp, j + offset, i) * G.Volume<F2>(k + c_half + ktemp, j + offset, i);
            //                }
            //                B_U(bl)(F3, 0, k + c_half, j_c, i) = (B_U(bl)(F3, 0, k, j_c, i) * G.Volume<F3>(k, j_c, i) 
            //                                                + B_U(bl)(F3, 0, k + 2 * c_half, j_c, i) * G.Volume<F3>(k + 2 * c_half, j_c, i) 
            //                                                + point_out * (c_right_v - c_left_v))
            //                                            / (2. * G.Volume<F3>(k + c_half, j_c, i));
            //                //if (m::abs(B_U(bl)(F3, 0, k + c_half, j_c, i)) < 1e-9) B_U(bl)(F3, 0, k + c_half, j_c, i) = 0.; // If very small, then just put it to zero TODO: check with Ben. Is this ok?
            //                if (k == 8 && i == 8 && j == j_f) printf("HYERIN: c_left_v %.8g c_right_v %.8g F3_1 %.8g vol_1 %.8g F3_2 %.8g vol %.8g F3 %.8g vol %.8g\n", c_left_v, c_right_v, B_U(bl)(F3, 0, k, j_c, i), G.Volume<F3>(k, j_c, i), B_U(bl)(F3, 0, k + 2 * c_half, j_c, i), G.Volume<F3>(k + 2 * c_half, j_c, i), B_U(bl)(F3, 0, k + c_half, j_c, i), G.Volume<F3>(k + c_half, j_c, i));
            //                // Then modify the internal faces for lower levels
            //                int lv_cell_len, k_i;
            //                if (c_half > 1) {
            //                    //for (int lv = 0; lv < c_half - 1; ++lv) {
            //                    for (int lv = c_half - 2; lv >= 0; --lv) {
            //                        lv_cell_len = m::pow(2, lv + 1); // cell length for each level of derefinement
            //                        for (int ntemp = 0; ntemp < (coarse_cell_len / lv_cell_len); ++ntemp) {
            //                            k_i = (2 * ntemp + 1) * lv_cell_len / 2; // internal face indices for each level
            //                            // take average on each end of the level-dependent coarse cell.
            //                            B_U(bl)(F3, 0, k + k_i, j_c, i) = (B_U(bl)(F3, 0, k + k_i - lv_cell_len / 2, j_c, i)
            //                                                                * G.Volume<F3>(k + k_i - lv_cell_len / 2, j_c, i)
            //                                                             + B_U(bl)(F3, 0, k + k_i + lv_cell_len / 2, j_c, i)
            //                                                                * G.Volume<F3>(k + k_i + lv_cell_len / 2, j_c, i)) 
            //                                                            / (2. * G.Volume<F3>(k + k_i, j_c, i));
            //                        }
            //                    }
            //                    if (k == 8 && i ==8 && j == j_f) printf("HYERIN: F3 internal %.3g\n", B_U(bl)(F3, 0, k + 1, j_c, i));
            //                }
            //            }
            //        }
            //    );
            //}
            
            // OLD
            //pmb0->par_for("B_CT_derefine_poles", block.s, block.e, bF1.ks, bF1.ke, jps, jpe, bF1.is, bF1.ie,
            //    KOKKOS_LAMBDA (const int &bl, const int &k, const int &j, const int &i) {
            //        const auto& G = B_U.GetCoords(bl);
            //        
            //        int coarse_cell_len = 2 * (((binner) ? j - j_f : j_f - j) + 1);
            //        int c_half = coarse_cell_len / 2; // half of the coarse cell's length
            //        if (k % coarse_cell_len == 0) {
            //            // F1: just average over the two fine cells
            //            //Real avg = (B_U(bl)(F1, 0, k, j_c, i) * G.Volume<F1>(k, j_c, i) + 
            //            //            B_U(bl)(F1, 0, k + 1, j_c, i) * G.Volume<F1>(k + 1, j_c, i)) / 2.;
            //            Real avg = 0.;
            //            auto j_c = j + (binner ? 0 : -1); // cell center
            //            for (int ktemp = 0; ktemp < coarse_cell_len; ++ktemp) 
            //                avg += B_U(bl)(F1, 0, k + ktemp, j_c, i) * G.Volume<F1>(k + ktemp, j_c, i);
            //            avg /= coarse_cell_len;
            //            for (int ktemp = 0; ktemp < coarse_cell_len; ++ktemp)
            //                B_U(bl)(F1, 0, k + ktemp, j_c, i) = avg / G.Volume<F1>(k + ktemp, j_c, i);
            //            //B_U(bl)(F1, 0, k, j_c, i) = avg / G.Volume<F1>(k, j_c, i);
            //            //B_U(bl)(F1, 0, k + 1, j_c, i) = avg / G.Volume<F1>(k + 1, j_c, i);
            //            //if (k == 4 && i == 15) printf("HYERIN: i,j,k = (%d, %d, %d) BF1 = %.3g, BF2 = %.3g\n", i,j,k, B_U(bl)(F1,0,k,j,i), B_U(bl)(F1,0,k+1,j,i));
            //            
            //            // F2: 
            //            //if (j == j_f) { //The two fine cells have 0 fluxes through the physical-ghost boundaries.
            //            //    for (int ktemp = 0; ktemp < coarse_cell_len; ++ktemp) 
            //            //        B_U(bl)(F2, 0, k + ktemp, j, i) = 0.;
            //            //} else { // average the cells
            //                avg = 0.;
            //                for (int ktemp = 0; ktemp < coarse_cell_len; ++ktemp)
            //                    avg += B_U(bl)(F2, 0, k + ktemp, j, i) * G.Volume<F2>(k + ktemp, j, i);
            //                avg /= coarse_cell_len;
            //                for (int ktemp = 0; ktemp < coarse_cell_len; ++ktemp)
            //                    B_U(bl)(F2, 0, k + ktemp, j, i) = avg / G.Volume<F2>(k + ktemp, j, i);
            //            //}

            //            // F3: The internal face between the two fine cells will take care of the divB=0. The two other faces of the coarse cell will remain unchanged.
            //            // First modify the very central internal face. In other words, deal with the highest level internal face first.
            //            // Sum of F2 fluxes in the left and right half of the coarse cell each. 
            //            Real c_left_v = 0., c_right_v = 0.;
            //            for (int ktemp = 0; ktemp < c_half; ++ktemp) {
            //                c_left_v += B_U(bl)(F2, 0, k + c_half - ktemp - 1, j + offset, i) * G.Volume<F2>(k + c_half - ktemp - 1, j + offset, i);
            //                c_right_v += B_U(bl)(F2, 0, k + c_half + ktemp, j + offset, i) * G.Volume<F2>(k + c_half + ktemp, j + offset, i);
            //            }
            //            B_U(bl)(F3, 0, k + c_half, j_c, i) = (B_U(bl)(F3, 0, k, j_c, i) * G.Volume<F3>(k, j_c, i) 
            //                                            + B_U(bl)(F3, 0, k + 2 * c_half, j_c, i) * G.Volume<F3>(k + 2 * c_half, j_c, i) 
            //                                            + point_out * (c_right_v - c_left_v))
            //                                        / (2. * G.Volume<F3>(k + c_half, j_c, i));
            //            // Then modify the internal faces for lower levels
            //            int lv_cell_len, k_i;
            //            if (c_half > 1) {
            //                //for (int lv = 0; lv < c_half - 1; ++lv) {
            //                for (int lv = c_half - 2; lv >= 0; --lv) {
            //                    lv_cell_len = m::pow(2, lv + 1); // cell length for each level of derefinement
            //                    for (int ntemp = 0; ntemp < (coarse_cell_len / lv_cell_len); ++ntemp) {
            //                        k_i = (2 * ntemp + 1) * lv_cell_len / 2; // internal face indices for each level
            //                        // take average on each end of the level-dependent coarse cell.
            //                        B_U(bl)(F3, 0, k + k_i, j_c, i) = (B_U(bl)(F3, 0, k + k_i - lv_cell_len / 2, j_c, i)
            //                                                            * G.Volume<F3>(k + k_i - lv_cell_len / 2, j_c, i)
            //                                                         + B_U(bl)(F3, 0, k + k_i + lv_cell_len / 2, j_c, i)
            //                                                            * G.Volume<F3>(k + k_i + lv_cell_len / 2, j_c, i)) 
            //                                                        / (2. * G.Volume<F3>(k + k_i, j_c, i));
            //                    }
            //                }
            //            }

            //            // average the fluid quantities TODO: separate this into a separate routine.
            //            avg = (rho_U(bl)(0, k, j_c, i) + rho_U(bl)(0, k + 1, j_c, i)) / 2.;
            //            rho_U(bl)(0, k, j_c, i) = avg;
            //            rho_U(bl)(0, k + 1, j_c, i) = avg;
            //            avg = (u_U(bl)(0, k, j_c, i) + u_U(bl)(0, k + 1, j_c, i)) / 2.;
            //            u_U(bl)(0, k, j_c, i) = avg;
            //            u_U(bl)(0, k + 1, j_c, i) = avg;
            //            VLOOP {
            //                avg = (uvec_U(bl)(v, k, j_c, i) + uvec_U(bl)(v, k + 1, j_c, i)) / 2.;
            //                uvec_U(bl)(v, k, j_c, i) = avg;
            //                uvec_U(bl)(v, k + 1, j_c, i) = avg;
            //            }
            //        }
            //    }
            //);
        }
    }
    return TaskStatus::complete;
}
