/*
 *  File: inverter.cpp
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
#include "inverter.hpp"

// inverter.hpp includes the template and instantiations in the correct order
#include "domain.hpp"
#include "floors_functions.hpp"
#include "flux.hpp"
#include "reductions.hpp"
#include "multizone.hpp"

int Inverter::CountPFlags(MeshData<Real> *md)
{
    return Reductions::CountFlags(md, "pflag", Inverter::status_names, IndexDomain::interior, false)[0];
}

std::shared_ptr<KHARMAPackage> Inverter::Initialize(ParameterInput *pin, std::shared_ptr<Packages_t>& packages)
{
    auto pkg = std::make_shared<KHARMAPackage>("Inverter");
    Params &params = pkg->AllParams();

    // Inversion scheme.  Could be separate packages but they do share a lot,
    // and could share more e.g. inline floor applications
    std::vector<std::string> allowed_inverter_names = {"none", "onedw", "kastaun", "mixed"};
    std::string inverter_name = pin->GetOrAddString("inverter", "type", "kastaun", allowed_inverter_names);
    bool use_kastaun = false;
    if (inverter_name == "onedw") {
        params.Add("inverter_type", Type::onedw);
    } else if (inverter_name == "kastaun") {
        params.Add("inverter_type", Type::kastaun);
        use_kastaun = true;
    } else if (inverter_name == "mixed") {
        params.Add("inverter_type", Type::mixed);
    } else if (inverter_name == "none") {
        params.Add("inverter_type", Type::none);
    }

    // Solver options
    // Any other Noble et al. implemented for fun should use lower tol/iter count, see Noble+06
    Real err_tol = pin->GetOrAddReal("inverter", "err_tol", (use_kastaun) ? 1e-12 : 1e-8);
    params.Add("err_tol", err_tol);
    int iter_max = pin->GetOrAddInteger("inverter", "iter_max", (use_kastaun) ? 25 : 8);
    params.Add("iter_max", iter_max);

    // Floor options
    // Use a custom block for inverter floors to allow customization.  Not sure anyone *wants* that but...
    if (!pin->DoesBlockExist("inverter_floors")) {
        params.Add("inverter_prescription", Floors::MakePrescription(pin, "floors"));
            if (pin->DoesBlockExist("floors_inner"))
                params.Add("inverter_prescription_inner", Floors::MakePrescriptionInner(pin, Floors::MakePrescription(pin, "floors"), "floors_inner"));
            else
                params.Add("inverter_prescription_inner", Floors::MakePrescriptionInner(pin, Floors::MakePrescription(pin, "floors"), "floors"));
    } else {
        params.Add("inverter_prescription", Floors::MakePrescription(pin, "inverter_floors"));
        params.Add("inverter_prescription_inner", Floors::MakePrescriptionInner(pin, Floors::MakePrescription(pin, "inverter_floors"), "inverter_floors"));
    }

    // Fixup options
    // Whether to apply Normal frame floors right after the inversion
    bool apply_floors_with_inversion = pin->GetOrAddBoolean("inverter", "apply_floors_with_inversion", use_kastaun);
    params.Add("apply_floors_with_inversion", apply_floors_with_inversion);
    // Tolerance in the momenta before a zone is considered failed and the unsolvable velocity zeroed out
    // Real bad_vel_tolerance = pin->GetOrAddReal("inverter", "bad_vel_tolerance", 1e-8);
    // params.Add("bad_vel_tolerance", bad_vel_tolerance);

    // Fix by averaging neighboring cells.  Enabled by default for 1Dw, but Kastaun failures are more dire
    bool fix_average_neighbors = pin->GetOrAddBoolean("inverter", "fix_average_neighbors", !use_kastaun);
    params.Add("fix_average_neighbors", fix_average_neighbors);
    // Fix by zeroing the velocity, for particularly nasty zones.  Applied immediately rather than in fixup
    // Generally you don't want this -- if you're eliminating velocity and therefore momentum, better to
    // eliminate inertia too by setting floor density/temp
    bool fix_zero_velocity = pin->GetOrAddBoolean("inverter", "fix_zero_velocity", false);
    params.Add("fix_zero_velocity", fix_zero_velocity);
    // Fix by replacing with floors, uvec=0. Backstop for states which are just impossible to use
    bool fix_atmosphere = pin->GetOrAddBoolean("inverter", "fix_atmosphere", true);
    params.Add("fix_atmosphere", fix_atmosphere);

    // Flag denoting UtoP inversion failures
    // Needs boundary sync if the fixup code will use neighbors, and if
    // we're syncing prims and fixing up after
    bool sync_prims = packages->Get("Driver")->Param<bool>("sync_prims");
    Metadata m;
    if (sync_prims && fix_average_neighbors) {
        m = Metadata({Metadata::Real, Metadata::Cell, Metadata::Derived, Metadata::OneCopy, Metadata::FillGhost});
    } else {
        m = Metadata({Metadata::Real, Metadata::Cell, Metadata::Derived, Metadata::OneCopy});
    }
    pkg->AddField("pflag", m);

    // When not using floors, we need to declare fflag for ourselves
    m = Metadata({Metadata::Real, Metadata::Cell, Metadata::Derived, Metadata::OneCopy, Metadata::Overridable});
    pkg->AddField("fflag", m);

    // This package may be loaded even when evolving implicitly, e.g. for FOFC
    // Only register our callbacks if they're needed for explicit evolution or a guess
    if (!pin->GetBoolean("GRMHD", "implicit") || pin->GetBoolean("emhd", "ideal_guess")) {
        // We exist basically to do this
        pkg->BlockUtoP = Inverter::BlockUtoP;
        // We want to run U->P on most boundaries when we're synchronizing conserved variables
        pkg->BoundaryUtoP = Inverter::BlockUtoP;
        // However, we apply domain boundaries to primitives.
        // Registering this additional function conveys that to the callers in `Packages` and `Boundaries`
        pkg->DomainBoundaryPtoU = Flux::BlockPtoUMHD;
    }

    pkg->PostStepDiagnosticsMesh = Inverter::PostStepDiagnostics;

    // List (vector) of HistoryOutputVars that will all be enrolled as output variables
    parthenon::HstVar_list hst_vars = {};
    // Count total floors as a history item
    hst_vars.emplace_back(parthenon::HistoryOutputVar(UserHistoryOperation::sum, CountPFlags, "PFlags"));
    // TODO entries for each individual flag?
    // add callbacks for HST output to the Params struct, identified by the `hist_param_key`
    pkg->AddParam<>(parthenon::hist_param_key, hst_vars);

    return pkg;
}

/**
 * Internal inversion fn, templated on inverter type.  Calls through to templated u_to_p
 * This is called with the correct template argument from BlockUtoP
 */
template<Inverter::Type inverter>
inline void BlockPerformInversion(MeshBlockData<Real> *rc, IndexDomain domain, bool coarse)
{
    auto pmb = rc->GetBlockPointer();

    PackIndexMap prims_map, cons_map;
    auto U = GRMHD::PackMHDCons(rc, cons_map);
    auto P = GRMHD::PackMHDPrims(rc, prims_map);
    const VarMap m_u(cons_map, true), m_p(prims_map, false);

    auto fflag = rc->PackVariables(std::vector<std::string>{"fflag"});
    auto pflag = rc->PackVariables(std::vector<std::string>{"pflag"});

    if (U.GetDim(4) == 0 || pflag.GetDim(4) == 0)
        return;

    const Real gam = pmb->packages.Get("GRMHD")->Param<Real>("gamma");

    auto &pars = pmb->packages.Get("Inverter")->AllParams();
    const Real err_tol = pars.Get<Real>("err_tol");
    const int iter_max = pars.Get<int>("iter_max");
    const bool apply_floors_with_inversion = pars.Get<bool>("apply_floors_with_inversion");
    //const Real bad_vel_tolerance = pars.Get<Real>("bad_vel_tolerance");
    const bool fix_zero_velocity = pars.Get<bool>("fix_zero_velocity");
    const Floors::Prescription inverter_floors       = pars.Get<Floors::Prescription>("inverter_prescription");
    const Floors::Prescription inverter_floors_inner = pars.Get<Floors::Prescription>("inverter_prescription_inner");
    const bool radius_dependent_floors = inverter_floors.radius_dependent_floors;

    const bool normal_frame_floors = (pmb->packages.AllPackages().count("Floors")) ?
                                        pmb->packages.Get("Floors")->Param<Floors::InjectionFrame>("frame") ==
                                            Floors::InjectionFrame::normal_kastaun :
                                        true;

    const auto& G = pmb->coords;
    const EMHD::EMHD_parameters& emhd_params = EMHD::GetEMHDParameters(pmb->packages);
    int active_iin = pmb->packages.Get("Multizone")->Param<int>("active_iin");

    // Get the primitives from our conserved versions
    // Notice by default, we recover variables for only the physical (interior or interior-ghost)
    // zones!  These are the only ones which are filled at our point in the step
    const IndexRange3 b = (domain == IndexDomain::entire)
                          ? KDomain::GetPhysicalRange(rc) : KDomain::GetRange(rc, domain, coarse);
    bool mixed_inverter = false;
    if constexpr (inverter == Inverter::Type::mixed) mixed_inverter = true;
    pmb->par_for("U_to_P", b.ks, b.ke, b.js, b.je, active_iin, b.ie,
        KOKKOS_LAMBDA (const int &k, const int &j, const int &i) {
            int pflagl;
            if (mixed_inverter) {
                if (G.r(k, j, i) < 1000) {
                    pflagl = Inverter::u_to_p<Inverter::Type::kastaun>(G, U, m_u, gam, k, j, i, P, m_p, Loci::center,
                                                        25, 1e-12);
                }
                else {
                    pflagl = Inverter::u_to_p<Inverter::Type::onedw>(G, U, m_u, gam, k, j, i, P, m_p, Loci::center,
                                                        iter_max, err_tol);
                }
            }
            else {
                pflagl = Inverter::u_to_p<inverter>(G, U, m_u, gam, k, j, i, P, m_p, Loci::center,
                                                        iter_max, err_tol);
            }
            // Apply floors immediately, attempting to correct bad inversions
            if (apply_floors_with_inversion) {
                Real rhoflr_max, uflr_max;

                int fflagl = 0; // Could read fflag here to preserve prior floors but there better not be any
                bool used_rho_to_slow = false;
                if (normal_frame_floors) {
                    fflagl |= Floors::determine_floors(G, P, m_p, gam, k, j, i, inverter_floors, inverter_floors_inner,
                        rhoflr_max, uflr_max);
                    // Add a floor to density which controls wayward velocities
                    if (fflagl && inverter_floors.use_rho_to_slow) {
                        // Calculate necessary rho
                        Real rho = std::max(P(m_p.RHO, k, j, i), rhoflr_max);
                        const Real u = std::max(P(m_p.UU, k, j, i), uflr_max);
                        Real h = 1. + gam * u / rho;
                        const Real alpha  = 1. / m::sqrt(-G.gcon(Loci::center, j, i, 0, 0));
                        const Real a_over_g = alpha / G.gdet(Loci::center, j, i);
                        FourVectors Dtmp;
                        GRMHD::calc_4vecs(G, P, m_p, k, j, i, Loci::center, Dtmp);
                        Real Scov[3] = {U(m_u.U1, k, j, i) * a_over_g - P(m_p.RHO, k, j, i) * Dtmp.ucov[1] * alpha,
                                        U(m_u.U2, k, j, i) * a_over_g - P(m_p.RHO, k, j, i) * Dtmp.ucov[2] * alpha,
                                        U(m_u.U3, k, j, i) * a_over_g - P(m_p.RHO, k, j, i) * Dtmp.ucov[3] * alpha};
                        Real Scon[3];
                        Real gupper[GR_DIM][GR_DIM], glower[GR_DIM][GR_DIM];
                        G.gcon(Loci::center, j, i, gupper);
                        G.gcov(Loci::center, j, i, glower);
                        Scon[0] = ((gupper[1][1] - gupper[0][1]*gupper[0][1]/gupper[0][0])*Scov[0] +
                                    (gupper[1][2] - gupper[0][1]*gupper[0][2]/gupper[0][0])*Scov[1] +
                                    (gupper[1][3] - gupper[0][1]*gupper[0][3]/gupper[0][0])*Scov[2]);

                        Scon[1] = ((gupper[2][1] - gupper[0][2]*gupper[0][1]/gupper[0][0])*Scov[0] +
                                    (gupper[2][2] - gupper[0][2]*gupper[0][2]/gupper[0][0])*Scov[1] +
                                    (gupper[2][3] - gupper[0][2]*gupper[0][3]/gupper[0][0])*Scov[2]);

                        Scon[2] = ((gupper[3][1] - gupper[0][3]*gupper[0][1]/gupper[0][0])*Scov[0] +
                                    (gupper[3][2] - gupper[0][3]*gupper[0][2]/gupper[0][0])*Scov[1] +
                                    (gupper[3][3] - gupper[0][3]*gupper[0][3]/gupper[0][0])*Scov[2]);
                        Real Ssq = 0.0;
                        SPACELOOP(ii) Ssq += Scon[ii] * Scov[ii];

                        Real gamma_max = inverter_floors.gamma_max;
                        const int radius_dependent_gamma_max = inverter_floors.radius_dependent_gamma_max;
                        if (radius_dependent_gamma_max > 0 && G.r(k, j, i) > 3) {
                            Real V02 = m::pow(inverter_floors.V0, 2.);
                            Real vchar2 = 1. / G.r(k, j, i) + 1. / Multizone::CalcRB(gam, inverter_floors.rs_bondi);
                            Real betagamma2_max = V02 * vchar2;
                            gamma_max = m::sqrt(betagamma2_max + 1.);
                        }
                        Real rhou0_old = U(m_u.RHO, k, j, i);
                        Real Ttt_old = U(m_u.UU, k, j, i); // - rhou0_old;
                        Real Ttr_old = U(m_u.U1, k, j, i);
                        Real Ttth_old = U(m_u.U2, k, j, i);
                        Real Ttphi_old = U(m_u.U2, k, j, i);
                        Real Ttrnet_old = U(m_u.U1, k, j, i) - P(m_p.RHO, k, j, i) * Dtmp.ucov[1] * G.gdet(Loci::center, j, i);
                        Real rho_old = P(m_p.RHO, k, j, i);
                        Real ucov1_old = Dtmp.ucov[1];
                        Real ucov2_old = Dtmp.ucov[2];
                        Real ucov3_old = Dtmp.ucov[3];

                        const Real rho_1 = m::sqrt(Ssq) / (gamma_max * m::sqrt(1. - 1. / (gamma_max * gamma_max)) * (h * gamma_max - alpha));
                        Real gamma = GRMHD::lorentz_calc(G, P, m_p, k, j, i, Loci::center);
                        if (gamma > gamma_max) {
                            fflagl |= Floors::FFlag::INVERTER_GAMMA;
                            used_rho_to_slow = true;
                            //rhoflr_max = rhoh_min -  gam * u;//rho * rhoh_min/rhoh;
                            //uflr_max = u; // * rhoh_min/rhoh;

                            Real Bvec[] = {0.0, 0.0, 0.0};
                            SPACELOOP(ii) Bvec[ii] = P(m_u.B1 + ii, k, j, i) * alpha;
                            Real BdotS = 0.;
                            SPACELOOP(ii) BdotS += Bvec[ii] * Scov[ii];
                            Real Bsq = 0.;
                            SPACELOOP2(ii, jj) Bsq += glower[ii + 1][jj + 1] * Bvec[ii] * Bvec[jj];
                            Real Sparsq = BdotS * BdotS / Bsq;
                            Real Sperpsq = Ssq - Sparsq;
                            Real Spar[3], Sperp[3];
                            SPACELOOP(ii) Spar[ii] = BdotS * Bvec[ii] / Bsq;
                            SPACELOOP(ii) Sperp[ii] = Scon[ii] - Spar[ii];

                    		// Equation for Lorentz factor W
                            auto frho = [&] (Real rhosol) {
                                const Real denom = rhosol * h * gamma_max * gamma_max - rhosol * alpha * gamma_max;
                                return Sparsq / SQR(denom)
                                        + Sperpsq / SQR(denom + Bsq)
                                        + 1./(gamma_max*gamma_max) - 1.;
                            };

                            Real zm = std::min(rho / 100., rho_1 / 100.);
                            Real zp = rho_1; //gamma_max;
                            Real z = 0.5*(zm + zp);
                            Real fm = frho(zm);
                            Real fp = frho(zp);

                            Real tol = 1e-8;
                            int iter;
                            for (iter = 0; iter < 30; ++iter) {
                                z =  (zm*fp - zp*fm)/(fp-fm);  // linear interpolation to point f(z)=0
                                Real f = frho(z);
                                // Quit if convergence reached
                                if ((m::abs(zm-zp) < tol) || (m::abs(f) < tol)) {
                                    break;
                                }
                                // root bracketed by [z,zp]
                                if (f*fp < 0.0) {
                                    zm = z;
                                    fm = f;
                                } else {  // root bracketed by [zm,z]
                                    zp = z;
                                    fp = f;
                                }
                            }
                            rho = z;
                            SPACELOOP(ii) P(m_p.U1+ii, k, j, i) = gamma_max * (Spar[ii] / (rho * h * gamma_max * gamma_max - rho * alpha * gamma_max) +
                                    Sperp[ii] / (rho * h * gamma_max * gamma_max + Bsq - rho * alpha * gamma_max));
                            GRMHD::calc_4vecs(G, P, m_p, k, j, i, Loci::center, Dtmp);
                            Real bsq = dot(Dtmp.bcon, Dtmp.bcov);
                            P(m_p.RHO, k, j, i) = rho;
                            P(m_p.UU, k, j, i) = (h - 1.) * rho / gam;
                            // P->U for any modified zones
                            Flux::p_to_u_mhd(G, P, m_p, emhd_params, gam, k, j, i, U, m_u, Loci::center);
                            Real rhou0_new = U(m_u.RHO, k, j, i);
                            Real Ttt_new = U(m_u.UU, k, j, i); // - rhou0_new;
                            Real Ttr_new = U(m_u.U1, k, j, i);
                            Real Ttrnet_new = U(m_u.U1, k, j, i) - P(m_p.RHO, k, j, i) * Dtmp.ucov[1] * G.gdet(Loci::center, j, i);

                            Real rhou0_diff = ((rhou0_new - rhou0_old) / rhou0_old);
                            Real Ttt_diff = ((Ttt_new - Ttt_old) / std::abs(Ttt_old));
                            Real Ttr_diff = ((Ttr_new - Ttr_old) / Ttr_old);
                            Real Ttrnet_diff = ((Ttrnet_new - Ttrnet_old) / Ttrnet_old);
                            Real mindiff = 1e-3;

                            if ((std::abs(Ttt_diff) > mindiff) || (std::abs(Ttr_diff) > mindiff))
                                printf("(i,j,k)=(%d,%d,%d), rhou0 = %.3g->%.3g (%.3g), Ttt = %.3g->%.3g (%.3g), Ttr = %.3g->%.3g (%.3g), Ttrnet = %.3g->%.3g (%.3g)\n", i,j,k, rhou0_old, rhou0_new, rhou0_diff, Ttt_old, Ttt_new, Ttt_diff, Ttr_old, Ttr_new, Ttr_diff, Ttrnet_old, Ttrnet_new, Ttrnet_diff);
                                //printf("(i,j,k)=(%d,%d,%d), rhoold=%.5g, ucov=(%.5g,%.5g,%.5g), alpha=%.5g, rhou0 = %.3g->%.3g (%.3g), Ttt = %.3g->%.3g (%.3g), Ttr = %.3g->%.3g (%.3g), Ttrnet = %.3g->%.3g (%.3g), Ttth=%.3g Ttphi%.3g\n", i,j,k, rho_old, ucov1_old, ucov2_old, ucov3_old, alpha, rhou0_old, rhou0_new, rhou0_diff, Ttt_old, Ttt_new, Ttt_diff, Ttr_old, Ttr_new, Ttr_diff, Ttrnet_old, Ttrnet_new, Ttrnet_diff, Ttth_old, Ttphi_old);
                        }
                    }
                } else {
                    // Bare minimum floors for numerics, before applying the rest in user-selected frame
                    rhoflr_max = inverter_floors.rho_min_const;
                    uflr_max = inverter_floors.u_min_const;
                    if (P(m_p.RHO, k, j, i) < rhoflr_max) fflagl |= Floors::FFlag::GEOM_RHO;
                    if (P(m_p.UU, k, j, i) < uflr_max) fflagl |= Floors::FFlag::GEOM_U;
                }
                if (fflagl && (!used_rho_to_slow)) {
                    // Apply floors to P -- this calls inversion again
                    pflagl = Floors::apply_floors<Floors::InjectionFrame::normal_kastaun>(G, P, m_p, gam, k, j, i,
                            rhoflr_max, uflr_max, U, m_u);
                    apply_ceilings(G, P, m_p, gam, k, j, i, inverter_floors, inverter_floors_inner, U, m_u);
                    // if still too low, apply bare minimum floors
                    rhoflr_max = inverter_floors.rho_min_const;
                    uflr_max = inverter_floors.u_min_const;
                    if (P(m_p.RHO, k, j, i) < rhoflr_max) {
                        fflagl |= Floors::FFlag::GEOM_RHO;
                        P(m_p.RHO, k, j, i) = rhoflr_max;
                    }
                    if (P(m_p.UU, k, j, i) < uflr_max) {
                        fflagl |= Floors::FFlag::GEOM_U;
                        P(m_p.UU, k, j, i) = uflr_max;
                    }
                    // P->U for any modified zones
                    Flux::p_to_u_mhd(G, P, m_p, emhd_params, gam, k, j, i, U, m_u, Loci::center);

                }
                
                // If we recovered the velocity, mark we used the gamma ceiling
                // It's just going to be the ceiling for no reason and mess everything up
                //const Real rho = P(m_p.RHO, k, j, i);
                //const Real u = P(m_p.UU, k, j, i);
                //const Real uvec[NVEC] = {P(m_p.U1, k, j, i), P(m_p.U2, k, j, i), P(m_p.U3, k, j, i)};
                //const Real B_P[NVEC] = {P(m_p.B1, k, j, i), P(m_p.B2, k, j, i), P(m_p.B3, k, j, i)};

                //Real rho_ut = 0., T[GR_DIM] = {0.};
                //GRMHD::p_to_u_mhd(G, rho, u, uvec, B_P, gam, k, j, i, rho_ut, T);
                //const Real& tol = bad_vel_tolerance;
                //if ((std::abs((T[1] - U(m_u.U1, k, j, i)) / U(m_u.U1, k, j, i)) > tol) ||
                //    (std::abs((T[2] - U(m_u.U2, k, j, i)) / U(m_u.U2, k, j, i)) > tol) ||
                //    (std::abs((T[3] - U(m_u.U3, k, j, i)) / U(m_u.U3, k, j, i)) > tol)) {
                //    P(m_p.U1, k, j, i) = 0.;
                //    P(m_p.U2, k, j, i) = 0.;
                //    P(m_p.U3, k, j, i) = 0.;
                //    fflagl |= Floors::FFlag::GAMMA;
                //}

                fflag(0, k, j, i) = fflagl;
            }
            pflag(0, k, j, i) = pflagl;
        }
    );
}

void Inverter::BlockUtoP(MeshBlockData<Real> *rc, IndexDomain domain, bool coarse)
{
    // This only chooses an implementation.  See BlockPerformInversion and implementations e.g. onedw.hpp
    auto& type = rc->GetBlockPointer()->packages.Get("Inverter")->Param<Type>("inverter_type");
    switch(type) {
    case Type::onedw:
        BlockPerformInversion<Type::onedw>(rc, domain, coarse);
        break;
    case Type::kastaun:
        BlockPerformInversion<Type::kastaun>(rc, domain, coarse);
        break;
    case Type::mixed:
        BlockPerformInversion<Type::mixed>(rc, domain, coarse);
        break;
    case Type::none:
        break;
    }
    // This is dangerous since there are many blocks/packs and we need one reduction. For later.
    //Reductions::StartFlagReduce(md, "pflag", Inverter::status_names, IndexDomain::interior, false, 1);
}

TaskStatus Inverter::PostStepDiagnostics(const SimTime& tm, MeshData<Real> *md)
{
    auto pmesh = md->GetMeshPointer();
    auto pmb0 = md->GetBlockData(0)->GetBlockPointer();
    // Options
    const auto& pars = pmesh->packages.Get("Globals")->AllParams();
    const int flag_verbose = pars.Get<int>("flag_verbose");

    // Debugging/diagnostic info about inversion flags
    // TODO grab the total and die on too many
    if (flag_verbose >= 1) {
        // TODO this should move into UtoP when everything goes MeshData
        Reductions::StartFlagReduce(md, "pflag", Inverter::status_names, IndexDomain::interior, false, 1);
        Reductions::CheckFlagReduceAndPrintHits(md, "pflag", Inverter::status_names, IndexDomain::interior, false, 1);

        // If we're the only floors, print those too
        if (!pmesh->packages.AllPackages().count("Floors")) {
            Reductions::StartFlagReduce(md, "fflag", Floors::FFlag::flag_names, IndexDomain::interior, true, 0);
            // Debugging/diagnostic info about floors
            Reductions::CheckFlagReduceAndPrintHits(md, "fflag", Floors::FFlag::flag_names, IndexDomain::interior, true, 0);
        }
    }

    return TaskStatus::complete;
}
