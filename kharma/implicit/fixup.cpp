/* 
 *  File: fixup.cpp
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

#include "implicit.hpp"

#include "floors.hpp"
#include "flux_functions.hpp"

TaskStatus Implicit::FixSolve(MeshBlockData<Real> *mbd) {

    Flag(mbd, "Fixing implicit solver failures");
    // Get MeshBlock pointer and obtain flag for primitives
    auto pmb = mbd->GetBlockPointer();

    // Get number of implicit variables
    PackIndexMap implicit_prims_map;
    auto implicit_vars = Implicit::GetOrderedNames(mbd, Metadata::GetUserFlag("Primitive"), true);
    auto& P            = mbd->PackVariables(implicit_vars, implicit_prims_map);
    const int nfvar    = P.GetDim(4);

    // Get grid object
    const auto& G = pmb->coords;

    GridScalar solve_fail = mbd->Get("solve_fail").data;
    GridScalar fflag      = mbd->Get("fflag").data;

    const Real gam    = pmb->packages.Get("GRMHD")->Param<Real>("gamma");
    // TODO flag_verbose here. Merge with other fixup into separate package or in GRMHD?
    // We'll want to try new in-depth fixes w/implicit as we go...
    const int verbose = pmb->packages.Get("Globals")->Param<int>("verbose");
    const Floors::Prescription floors(pmb->packages.Get("Floors")->AllParams());

    // Boundaries were synced just before the call to this function (cf. imex_driver.cpp). 
    // Which means unsuccessful values were copied to ghost zones. Therefore, we need to loop over entire domain.
    const IndexRange ib = mbd->GetBoundsI(IndexDomain::entire);
    const IndexRange jb = mbd->GetBoundsJ(IndexDomain::entire);
    const IndexRange kb = mbd->GetBoundsK(IndexDomain::entire);

    auto bounds  = pmb->cellbounds;
    const int n1 = bounds.ncellsi(IndexDomain::entire);
    const int n2 = bounds.ncellsj(IndexDomain::entire);
    const int n3 = bounds.ncellsk(IndexDomain::entire);

    const IndexRange ib_b = mbd->GetBoundsI(IndexDomain::interior);
    const IndexRange jb_b = mbd->GetBoundsJ(IndexDomain::interior);
    const IndexRange kb_b = mbd->GetBoundsK(IndexDomain::interior);

    // TODO don't allocate here
    ParArrayND<Real> sum("sum_good_neighbors", nfvar, n3+1, n2+1, n1+1);
    ParArrayND<Real> sum_x("sum_all_neighbors", nfvar, n3+1, n2+1, n1+1);

    pmb->par_for("fix_solver_failures", kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
        KOKKOS_LAMBDA (const int& k, const int& j, const int& i) {
            FLOOP {
                sum(ip, k, j, i)   = 0.;
                sum_x(ip, k, j, i) = 0.;
            }
            // Fix only bad zones
            if ((solve_fail(k, j, i)) == SolverStatus::fail) {
                //printf("Fixing zone %d %d %d!\n", i, j, k);
                double wsum = 0., wsum_x = 0.;
                // double sum[nfvar] = {0.}, sum_x[nfvar] = {0.};
                // For all neighboring cells...
                for (int n = -1; n <= 1; n++) {
                    for (int m = -1; m <= 1; m++) {
                        for (int l = -1; l <= 1; l++) {
                            int ii = i + l, jj = j + m, kk = k + n;
                            // If we haven't overstepped array bounds...
                            if (inside(kk, jj, ii, kb, jb, ib)) {
                                // Weight by distance
                                // TODO abs(l) == l*l always?
                                double w = 1./(m::abs(l) + m::abs(m) + m::abs(n) + 1);

                                // Count only the good cells, if we can
                                if ((solve_fail(kk, jj, ii)) != SolverStatus::fail) {
                                    // Weight by distance.  Note interpolated "fixed" cells stay flagged
                                    wsum += w;
                                    FLOOP sum(ip, k, j, i) += w * P(ip, kk, jj, ii);
                                }
                                // Just in case, keep a sum of even the bad ones
                                wsum_x += w;
                                FLOOP sum_x(ip, k, j, i) += w * P(ip, kk, jj, ii);
                            }
                        }
                    }
                }

                if(wsum < 1.e-10) {
                    // TODO probably should crash here. Or average anyway?
#ifndef KOKKOS_ENABLE_SYCL
                    if (verbose >= 1 && inside(k, j, i, kb_b, jb_b, ib_b)) // If an interior zone...
                        printf("No neighbors were available at %d %d %d!\n", i, j, k);
#endif // TODO SYCL has cout
                } else {
                    FLOOP P(ip, k, j, i) = sum(ip, k, j, i)/wsum;
                }
            }
        }
    );

    // Since floors were applied earlier, we assume the zones obtained by averaging the neighbors also respect the floors.
    // Compute new conserved variables
    PackIndexMap prims_map, cons_map;
    auto& P_all = mbd->PackVariables(std::vector<MetadataFlag>{Metadata::GetUserFlag("Primitive")}, prims_map);
    auto& U_all = mbd->PackVariables(std::vector<MetadataFlag>{Metadata::Conserved}, cons_map);
    const VarMap m_u(cons_map, true), m_p(prims_map, false);
    // Get new sizes
    const int nvar = P_all.GetDim(4);

    // Need emhd_params object
    EMHD_parameters emhd_params = EMHD::GetEMHDParameters(pmb->packages);

    pmb->par_for("fix_solver_failures_PtoU", kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
        KOKKOS_LAMBDA (const int& k, const int& j, const int& i) {
            if (( solve_fail(k, j, i)) == SolverStatus::fail)
                Flux::p_to_u(G, P_all, m_p, emhd_params, gam, k, j, i, U_all, m_u);
        }
    );

    Flag(mbd, "Fixed solver failures");
    return TaskStatus::complete;

}
