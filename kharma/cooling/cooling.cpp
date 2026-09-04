/* 
 *  File: cooling.cpp
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
#include "cooling.hpp"
#include "grmhd.hpp"

std::shared_ptr<KHARMAPackage> Cooling::Initialize(ParameterInput *pin, std::shared_ptr<Packages_t>& packages)
{
    auto pkg = std::make_shared<KHARMAPackage>("Cooling");
    Params &params = pkg->AllParams();

    // Cooling term
    std::vector<std::string> cooling_options = {"noble"};
    std::string cooling_type = pin->GetOrAddString("cooling", "type", "noble", cooling_options);
    params.Add("cooling_type", cooling_type);
    Real target_entropy = pin->GetOrAddReal("cooling", "target_entropy", 0.01);
    params.Add("target_entropy", target_entropy);

    pkg->AddSource = Cooling::AddSource;

    return pkg;
}

TaskStatus Cooling::AddSource(MeshData<Real> *md, MeshData<Real> *mdudt, IndexDomain domain)
{
    // Pointers
    auto pmesh = mdudt->GetMeshPointer();
    auto pmb0 = mdudt->GetBlockData(0)->GetBlockPointer();
    // Options
    const auto& gpars = pmb0->packages.Get("GRMHD")->AllParams();
    const auto& pars = pmb0->packages.Get("Cooling")->AllParams();

    // Pack variables
    PackIndexMap prims_map, cons_map;
    auto dUdt = mdudt->PackVariables(std::vector<MetadataFlag>{Metadata::Conserved}, cons_map);
    auto P    = md->PackVariables(std::vector<MetadataFlag>{Metadata::GetUserFlag("Primitive")}, prims_map);
    const VarMap m_u(cons_map, true), m_p(prims_map, false);
    // Get sizes
    const IndexRange ib = mdudt->GetBoundsI(IndexDomain::interior);
    const IndexRange jb = mdudt->GetBoundsJ(IndexDomain::interior);
    const IndexRange kb = mdudt->GetBoundsK(IndexDomain::interior);
    const IndexRange block = IndexRange{0, dUdt.GetDim(5) - 1};
    const Real gam = gpars.Get<Real>("gamma");
    const Real Sstar = pars.Get<Real>("target_entropy");
    

    pmb0->par_for("add_cooling", block.s, block.e, kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
        KOKKOS_LAMBDA (const int& b, const int &k, const int &j, const int &i) {
            const auto& G = dUdt.GetCoords(b);
            // Need coordinates to evaluate particle addtn rate
            // Note that makes the wind spherical-only, TODO ensure this
            GReal Xembed[GR_DIM];
            G.coord_embed(k, j, i, Loci::center, Xembed);
            GReal r = Xembed[1], th = Xembed[2];
            GReal a = G.coords.get_a();

            // Following eq.7 of Avara+24
            GReal rho = P(b, m_p.RHO, k, j, i);
            GReal u = P(b, m_p.UU, k, j, i);
            GReal Omega = 1. / (a + SQR(r * r * r));
            GReal torb = 2. * M_PI / Omega;
            GReal S = (gam - 1.) * u / m::pow(rho, gam); // entropy
            GReal L = u * SQR(S / Sstar - 1.) / torb;
            
            // evaluate Be
            FourVectors D;
            GRMHD::calc_4vecs(G, P(b), m_p, k, j, i, Loci::center, D);
            GReal bsq = dot(D.bcon, D.bcov);
            GReal Be = -(1. + (gam * u + bsq) / rho) - 1.;// Bernoulli parameter (Penna+13)

            if ((S <= Sstar) && (Be > 0)) L = 0.;

            Real new_du[GR_DIM] = {0};
            for (int lam = 0; lam < GR_DIM; ++lam)
                new_du[lam] += - G.gdet(Loci::center, j, i) * D.ucov[lam] * L;
            
            dUdt(b, m_u.UU, k, j, i)           += new_du[0];
            VLOOP dUdt(b, m_u.U1 + v, k, j, i) += new_du[1 + v];
        }
    );

    return TaskStatus::complete;
}
