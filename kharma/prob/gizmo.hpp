/* 
 *  File: bondi.hpp
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
#pragma once

#include "decs.hpp"

#include "bondi.hpp"
#include "gr_coordinates.hpp"
#include "flux_functions.hpp"
#include "grmhd_functions.hpp"
#include "pack.hpp"
#include "coordinate_utils.hpp"
#include "types.hpp"

#include <parthenon/parthenon.hpp>

/**
 * Initialize a Bondi problem over the domain
 */
TaskStatus InitializeGIZMO(std::shared_ptr<MeshBlockData<Real>>& rc, ParameterInput *pin);

/**
 * Set all values on a given domain to the Bondi inflow analytic steady-state solution
 * 
 * Used for initialization and boundary conditions
 */
TaskStatus SetGIZMO(std::shared_ptr<MeshBlockData<Real>>& rc, IndexDomain domain, bool coarse=false);

KOKKOS_INLINE_FUNCTION void XtoindexGIZMO(const GReal XG[GR_DIM],
                                    const GridScalar& rarr, const int length, int& i, GReal& del)
{
    Real dx2, dx2_min;
    dx2_min = m::pow(XG[1]-rarr(0),2); //100000.; //arbitrarily large number

    i = 0; // initialize

    for (int itemp = 0; itemp < length; itemp++) {
        if (rarr(itemp) < XG[1]) { // only look for smaller side
            dx2 = m::pow(XG[1] - rarr(itemp), 2);

            // simplest interpolation (Hyerin 07/26/22)
            if (dx2 < dx2_min){
                dx2_min = dx2;
                i = itemp;
            }
        }
    }
    
    // interpolation (11/14/2022) TODO: write a case where indices hit the boundaries of the data file
    del = (XG[1]-rarr(i))/(rarr(i+1)-rarr(i));

    if (m::abs(dx2_min/m::pow(XG[1],2))>1.e-8) printf("XtoindexGizmo: dx2 pretty large = %g at r= %g \n",dx2_min, XG[1]);
}
/**
 * Get the GIZMO output values at a particular zone
 * Note this assumes that there are ghost zones!
 */
KOKKOS_INLINE_FUNCTION void get_prim_gizmo_shell(const GRCoordinates& G, const VariablePack<Real>& P, const VarMap& m_p, 
                                           const Real& gam, const Real rin_init, const Real rs, const Real ur_frac, const Real uphi,
                                           const GridScalar& rarr, const GridScalar& rhoarr, const GridScalar& Tarr, const GridScalar& vrarr, const int length,
                                           const int& k, const int& j, const int& i)
{
    // Solution constants for velocity prescriptions
    // Ideally these could be cached but preformance isn't an issue here
    Real mdot = 1.; // mdot defined arbitrarily

    GReal Xnative[GR_DIM], Xembed[GR_DIM];
    G.coord(k, j, i, Loci::center, Xnative);
    G.coord_embed(k, j, i, Loci::center, Xembed);
    GReal r = Xembed[1];

    // Get GIZMO or Bondi data
    Real rho, u, ur;
    if (r < rin_init * 0.9){
        // Bondi solution
        get_bondi_soln(r, rs, mdot, gam, rho, u, ur);
    } else {
        // linear interpolation
        int itemp; GReal del;
        XtoindexGIZMO(Xembed, rarr, length, itemp, del);
        if (del < 0 ) { // when r is smaller than GIZMO's range
            del = 0; // just copy over the smallest r values
        }
        rho = rhoarr(itemp) * (1.-del) + rhoarr(itemp+1) * del;
        u = rho * (Tarr(itemp) * (1.-del) + Tarr(itemp+1) * del) / (gam - 1.);
        ur = 0.;
    }
    Real ucon_bl[GR_DIM] = {0., ur * ur_frac, 0., uphi * m::pow(r, -3. / 2.)};

    // Set u^t and transform to native coordinates
    GReal ucon_native[GR_DIM];
    G.coords.bl_fourvel_to_native(Xnative, ucon_bl, ucon_native);

    // Convert native 4-vector to primitive u-twiddle, see Gammie '04
    Real gcon[GR_DIM][GR_DIM], u_prim[NVEC];
    G.gcon(Loci::center, j, i, gcon);
    fourvel_to_prim(gcon, ucon_native, u_prim);

    P(m_p.RHO, k, j, i) = rho;
    P(m_p.UU, k, j, i) = u;
    P(m_p.U1, k, j, i) = u_prim[0];
    P(m_p.U2, k, j, i) = u_prim[1];
    P(m_p.U3, k, j, i) = u_prim[2];
}

KOKKOS_INLINE_FUNCTION Real frac_diff(const GReal XG[GR_DIM], const GReal x1, const GReal x2, const GReal x3) {
    //Real dx2 = m::pow((x1 - XG[1]) / XG[1], 2.) + m::pow((x2 - XG[2]) / M_PI, 2.) + m::pow((x3 - XG[3]) / (2. * M_PI), 2.);
    Real phi_temp;
    phi_temp = (XG[3] > 2 * M_PI)? XG[3] - 2 * M_PI : XG[3];
    phi_temp = (phi_temp < 0)? phi_temp + 2 * M_PI : phi_temp;
    Real dx2 = m::pow((x1 - XG[1]), 2.) + m::pow((x2 - XG[2]) / M_PI, 2.) + m::pow((x3 - phi_temp) / (2. * M_PI), 2.);
    return dx2;
}

KOKKOS_INLINE_FUNCTION void XtoindexGIZMO3D(const GReal XG[GR_DIM],
                                    const GridVector& coordarr, const hsize_t length, int& i, GReal& del)
{
    Real dx2, dx2_min;
    dx2_min = frac_diff(XG, coordarr(0,0), coordarr(0,1), coordarr(0,2));

    i = 0; // initialize

    for (int itemp = 0; itemp < length; itemp++) {
        if (m::abs(dx2_min) > 1.0e-8) { // continue if the difference is large
            dx2 = frac_diff(XG, coordarr(itemp, 0), coordarr(itemp, 1), coordarr(itemp, 2));

            // simplest interpolation (Hyerin 07/26/22)
            if (dx2 < dx2_min){
                dx2_min = dx2;
                i = itemp;
            }
        }
    }
    
    // No interpolation! Warn if the data points are not exactly on top of each other
    if ((m::abs(dx2_min) > 1.e-8) && (XG[2] > 0) && (XG[2] < M_PI)) printf("XtoindexGIZMO3D: dx2 frac diff large = %g at (r,th,phi)=(%lf %lf %lf) fitted=(%lf %lf %lf) i = %d \n",m::sqrt(dx2_min), XG[1], XG[2], XG[3], coordarr(i,0),coordarr(i,1),coordarr(i,2), i);
}

/**
 * Get the GIZMO output values at a particular zone for 3D GIZMO data
 * Note this assumes that there are ghost zones!
 */
KOKKOS_INLINE_FUNCTION void get_prim_gizmo_shell_3d(const GRCoordinates& G, const VariablePack<Real>& P, const VarMap& m_p,
                                           const Real& gam, const Real rin_init, const Real rs, const Real ur_frac, const Real uphi_frac,
                                           const GridVector& coordarr, const GridScalar& rhoarr, const GridScalar& Tarr, const GridVector& varr, const hsize_t length,
                                           const int& k, const int& j, const int& i)
{
    Real mdot = 1.; // mdot defined arbitrarily
    Real n = 1. / (gam - 1.);
    GReal Xnative[GR_DIM], Xembed[GR_DIM];//, Xembed_corner[GR_DIM];
    G.coord(k, j, i, Loci::center, Xnative);
    G.coord_embed(k, j, i, Loci::center, Xembed);
    //G.coord_embed(k, j, i, Loci::corner, Xembed_corner); // TODO: get cell centered values from KungYi
    GReal r = Xembed[1];
    GReal th = Xembed[2];

    Real rho, u, ur, uth, uphi;
    GReal del;
    if (r < rin_init * 0.9){
        // Bondi solution
        get_bondi_soln(r, rs, mdot, gam, rho, u, ur);
        ur = ur * ur_frac;
        uth = 0.0;
        uphi = 0.0;
    } else {
        int itemp;
        XtoindexGIZMO3D(Xembed, coordarr, length, itemp, del);
        // DO NOT INTERPOLATE, it is assumed GIZMO data is right on the grid
        rho = rhoarr(itemp);
        u = rho * (Tarr(itemp)) * n;
        ur = varr(itemp, 0);
        uth = varr(itemp, 1) / r;
        uphi = varr(itemp, 2) / (r * m::sin(th));
        if (r < 2e6) uphi = uphi_frac * m::pow(r, -3./2.); // override because below this value we don't trust GIZMO
    }
    
    Real ucon_bl[GR_DIM] = {0., ur, uth, uphi};

    // Set u^t and transform to native coordinates
    GReal ucon_native[GR_DIM];
    G.coords.bl_fourvel_to_native(Xnative, ucon_bl, ucon_native);

    // Convert native 4-vector to primitive u-twiddle, see Gammie '04
    Real gcon[GR_DIM][GR_DIM], u_prim[NVEC];
    G.gcon(Loci::center, j, i, gcon);
    fourvel_to_prim(gcon, ucon_native, u_prim);

    P(m_p.RHO, k, j, i) = rho;
    P(m_p.UU, k, j, i) = u;
    P(m_p.U1, k, j, i) = u_prim[0];
    P(m_p.U2, k, j, i) = u_prim[1];
    P(m_p.U3, k, j, i) = u_prim[2];

}


KOKKOS_INLINE_FUNCTION void get_B_gizmo_shell_3d(const GRCoordinates& G, const VariablePack<Real>& P,
                   const Real rin_init, 
                   const GridVector& coordarr, const GridVector& varr, const GridVector& B, const GridVector& B_save, const hsize_t length,
                   const int& k, const int& j, const int& i)
{
    //Real B_cons[NVEC];
    GReal Xnative[GR_DIM], Xembed[GR_DIM];
    G.coord(k, j, i, Loci::center, Xnative);
    G.coord_embed(k, j, i, Loci::center, Xembed);
    GReal r = Xembed[1];
    GReal th = Xembed[2];

    int itemp;
    GReal del;
    
    if (r < rin_init * 0.9){
        // don't do anything and wait for seed_B
    } else {
        XtoindexGIZMO3D(Xembed, coordarr, length, itemp, del);
        
        // First get velocities in embedded coordinates, taking Newtonian approximation
        Real ucon_embed[GR_DIM], ucov_embed[GR_DIM], ucon_native[GR_DIM]; //, gcov_ks[GR_DIM][GR_DIM]
        ucon_embed[1] = varr(itemp, 0);
        ucon_embed[2] = varr(itemp, 1) / r;
        ucon_embed[3] = varr(itemp, 2) / (r * m::sin(th));
        G.coords.ks_lower_fourvel(Xembed, ucon_embed, ucov_embed);
        //ks.gcov_embed(Xembed, gcov_ks); // gcov_embed
        //set_ut(gcov_ks, ucon_embed); // Set u^t to make a 4-vector

        //DLOOP1 ucov_embed[mu] = 0.;
        //DLOOP2 ucov_embed[mu] += gcov_ks[mu][nu] * ucon_embed[nu]; // lower
        G.coords.con_vec_to_native(Xnative, ucon_embed, ucon_native); // prepare ucon_native for later use

        // make magnetic fields into a four-vector
        Real B_embed[NVEC], B_native[NVEC], bcon_embed[GR_DIM], bcov_ks[GR_DIM], bcon_native[GR_DIM];
        B_embed[0] = B(itemp, 0);
        B_embed[1] = B(itemp, 1) / r;
        B_embed[2] = B(itemp, 2) / (r * m::sin(th));
        bcon_embed[0] = 0;
        VLOOP bcon_embed[0] += B_embed[v] * ucov_embed[v+1];
        VLOOP bcon_embed[v+1] = (B_embed[v] + bcon_embed[0] * ucon_embed[v+1]) / ucon_embed[0];

        // convert into native coordinates
        G.coords.con_vec_to_native(Xnative, bcon_embed, bcon_native);

        // convert the 4-vector into 3-vector
        VLOOP B_native[v] = bcon_native[v+1] * ucon_native[0] - bcon_native[0] * ucon_native[v+1];
        
        VLOOP B_save(v, k, j, i) = B_native[v] * G.gdet(Loci::center, j, i); // convert to B_U
    }

}
