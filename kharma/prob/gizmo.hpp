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
#include "hdf5_utils.h"

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

KOKKOS_INLINE_FUNCTION Real frac_diff(const GReal XG[GR_DIM], const GReal x1, const GReal x2, const GReal x3) {
    Real phi_temp;
    phi_temp = (XG[3] > 2 * M_PI)? XG[3] - 2 * M_PI : XG[3];
    phi_temp = (phi_temp < 0)? phi_temp + 2 * M_PI : phi_temp;
    Real dx2 = m::pow((x1 - XG[1]), 2.) + m::pow((x2 - XG[2]) / M_PI, 2.) + m::pow((x3 - phi_temp) / (2. * M_PI), 2.);
    return dx2;
}

KOKKOS_INLINE_FUNCTION void XtoindexGizmo3D(const GReal XG[GR_DIM],
                                    const GridVector& coordarr, const hsize_t length, int& i, GReal& del)
{
    Real dx2, dx2_min;
    dx2_min = frac_diff(XG, coordarr(0,0), coordarr(0,1), coordarr(0,2));

    i = 0; // initialize

    for (int itemp = 0; itemp < length; itemp++) {
        if (m::abs(dx2_min) > 1.0e-8) { // continue if the difference is large
            dx2 = frac_diff(XG, coordarr(itemp, 0), coordarr(itemp, 1), coordarr(itemp, 2));


            // simplest interpolation (Hyerin 07/26/22)
            if (dx2<dx2_min){
                dx2_min=dx2;
                i = itemp;
            }
        }
    }
    
    // No interpolation! Warn if the data points are not exactly on top of each other
    if ((m::abs(dx2_min)>1.e-8) && (XG[2] > 0) && (XG[2] < M_PI)) printf("XtoindexGizmo3D: dx2 frac diff large = %g at (r,th,phi)=(%lf %lf %lf) fitted=(%lf %lf %lf) i = %d \n",m::sqrt(dx2_min), XG[1], XG[2], XG[3], coordarr(i,0),coordarr(i,1),coordarr(i,2), i);
}
/**
 * Get the GIZMO output values at a particular zone
 * Note this assumes that there are ghost zones!
 */
KOKKOS_INLINE_FUNCTION void get_prim_gizmo_shell(const GRCoordinates& G, const CoordinateEmbedding& coords, const VariablePack<Real>& P, const VarMap& m_p,
                                           const Real& gam,
                                           const Real rin_init, const Real rs, Real vacuum_rho, Real vacuum_u_over_rho,
                                           const GridScalar& rarr, const GridScalar& rhoarr, const GridScalar& Tarr, const GridScalar& vrarr, const int length,
                                           const int& k, const int& j, const int& i)
{
    // Solution constants for velocity prescriptions
    // Ideally these could be cached but preformance isn't an issue here
    Real mdot = 1.; // mdot defined arbitrarily
    //Real rs = 1./sqrt(T); //1000.;

    GReal Xnative[GR_DIM], Xembed[GR_DIM];
    G.coord(k, j, i, Loci::center, Xnative);
    G.coord_embed(k, j, i, Loci::center, Xembed);
    GReal r = Xembed[1];

    // Get GIZMO or vacuum/Bondi data
    Real rho, u, ur;
    if (r < rin_init * 0.9){
        // Vacuum values for interior
        rho = vacuum_rho;
        u = vacuum_rho * vacuum_u_over_rho;
        // Radial velocity from Bondi solution
        Real rho_tmp, u_tmp;
        get_bondi_soln(r, rs, mdot, gam, rho_tmp, u_tmp, ur);
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
    Real ucon_bl[GR_DIM] = {0., ur, 0., 0.};

    // Set u^t and transform to native coordinates
    GReal ucon_native[GR_DIM];
    G.coords.bl_fourvel_to_native(Xnative, ucon_bl, ucon_native);

    // Convert native 4-vector to primitive u-twiddle, see Gammie '04
    Real gcon[GR_DIM][GR_DIM], u_prim[NVEC];
    G.gcon(Loci::center, j, i, gcon);
    fourvel_to_prim(gcon, ucon_native, u_prim);

    if(!isnan(rho))P(m_p.RHO, k, j, i) = rho;
    if(!isnan(u)) P(m_p.UU, k, j, i) = u;
    if(!isnan(u_prim[0])) P(m_p.U1, k, j, i) = u_prim[0];
    if(!isnan(u_prim[1])) P(m_p.U2, k, j, i) = u_prim[1];
    if(!isnan(u_prim[2])) P(m_p.U3, k, j, i) = u_prim[2];
}
/**
 * Get the GIZMO output values at a particular zone for 3D GIZMO data
 * Note this assumes that there are ghost zones!
 * TODO: Hyerin: maybe combine with get_prim_bondi and get_prim_gizmo_shell
 */
KOKKOS_INLINE_FUNCTION void get_prim_gizmo_shell_3d(const GRCoordinates& G, const CoordinateEmbedding& coords, const VariablePack<Real>& P, const VarMap& m_p,
                                           const Real& gam, const Real mdot, const Real rs, 
                                           const Real r_shell, const Real ur_frac, const Real uphi_frac, Real vacuum_logrho, Real vacuum_log_u_over_rho,
                                           const GridVector& coordarr, const GridScalar& rhoarr, const GridScalar& Tarr, const GridVector& varr, const hsize_t length,
                                           const int& k, const int& j, const int& i)
{
    GReal Xnative[GR_DIM], Xembed[GR_DIM];//, Xembed_corner[GR_DIM];
    G.coord(k, j, i, Loci::center, Xnative);
    G.coord_embed(k, j, i, Loci::center, Xembed);
    //G.coord_embed(k, j, i, Loci::corner, Xembed_corner); // TODO: get cell centered values from KungYi
    GReal r = Xembed[1];
    GReal th = Xembed[2];
    Real rho, u;
    Real u_prim[NVEC];

    if (r<r_shell){
        // leave it to get_prim_bondi
        get_prim_bondi(G, false, rs, mdot, gam, ur_frac, uphi_frac, 0., false, false, rho, u, u_prim, k, j, i);
        // override rho, u
        rho = pow(10.,vacuum_logrho); // pow(10.,-4.);
        u = rho * pow(10.,vacuum_log_u_over_rho);
        //ucon_bl[1] = ur_frac * ur;
        //ucon_bl[3] = uphi_frac * m::pow(r, -3./2.) * m::sin(th);
    } else {
        // Solution constants for velocity prescriptions
        // Ideally these could be cached but preformance isn't an issue here
        Real mdot = 1.; // mdot and rs defined arbitrarily
        Real n = 1. / (gam - 1.);
        Real uc = sqrt(mdot / (2. * rs));
        Real Vc = -sqrt(pow(uc, 2) / (1. - 3. * pow(uc, 2)));
        Real Tc = -n * pow(Vc, 2) / ((n + 1.) * (n * pow(Vc, 2) - 1.));
        Real C1 = uc * pow(rs, 2) * pow(Tc, n);
        Real C2 = pow(1. + (1. + n) * Tc, 2) * (1. - 2. * mdot / rs + pow(C1, 2) / (pow(rs, 4) * pow(Tc, 2 * n)));

        //Real smallrho=pow(10.,vacuum_logrho); // pow(10.,-4.);
        //Real smallu = smallrho*pow(10.,vacuum_log_u_over_rho);
        Real T, ur, uth, uphi;
        int itemp;
        GReal del;


        // also careful when initializing near the horizon
        if (r < 2) get_prim_bondi(G, false, rs, mdot, gam, ur_frac, uphi_frac, 0., false, false, rho, u, u_prim, k, j, i);
        else {
            T = get_T(r, C1, C2, n, rs);
            ur = -C1 / (pow(T, n) * pow(r, 2));
            Real ucon_bl[GR_DIM] = {0, 0, 0, 0};


            XtoindexGizmo3D(Xembed, coordarr, length, itemp, del);
            // DO NOT INTERPOLATE, it is assumed GIZMO data is right on the grid
            rho = rhoarr(itemp);
            u = rho * (Tarr(itemp)) * n;
            ur = varr(itemp, 0);
            uth = varr(itemp, 1) / r;
            uphi = varr(itemp, 2) / (r * m::sin(th));
            // Newtonian limit
            ucon_bl[1] = ur;
            ucon_bl[2] = uth;
            ucon_bl[3] = uphi;
            //if (r < 2e6) ucon_bl[3] = uphi_frac * m::pow(r, -3./2.) * m::sin(th); // override because below this value we don't trust GIZMO

            // Get the native-coordinate 4-vector corresponding to ur
            Real ucon_mks[GR_DIM];
            G.coords.bl_fourvel_to_native(Xnative, ucon_bl, ucon_mks);

            // Convert native 4-vector to primitive u-twiddle, see Gammie '04
            Real gcon[GR_DIM][GR_DIM];
            G.gcon(Loci::center, j, i, gcon);
            fourvel_to_prim(gcon, ucon_mks, u_prim);

        }
    }
    if(!isnan(rho)) P(m_p.RHO, k, j, i) = rho;
    if(!isnan(u)) P(m_p.UU, k, j, i) = u;
    if(!isnan(u_prim[0])) P(m_p.U1, k, j, i) = u_prim[0];
    if(!isnan(u_prim[1])) P(m_p.U2, k, j, i) = u_prim[1];
    if(!isnan(u_prim[2])) P(m_p.U3, k, j, i) = u_prim[2];

}

KOKKOS_INLINE_FUNCTION void get_B_gizmo_shell_3d(const GRCoordinates& G, const CoordinateEmbedding& coords, const VariablePack<Real>& P,
                   const Real r_shell, 
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
    
    if (r<r_shell){
        // don't do anything and wait for seed_B
    } else {
        XtoindexGizmo3D(Xembed, coordarr, length, itemp, del);
        
        // First get velocities in embedded coordinates, taking Newtonian approximation
        Real ucon_embed[GR_DIM], ucov_embed[GR_DIM], ucon_native[GR_DIM];
        ucon_embed[1] = varr(itemp, 0);
        ucon_embed[2] = varr(itemp, 1) / r;
        ucon_embed[3] = varr(itemp, 2) / (r * m::sin(th));
        Real gcov_bl[GR_DIM][GR_DIM];
        SphBLCoords(G.coords.get_a()).gcov_embed(Xembed, gcov_bl);
        set_ut(gcov_bl, ucon_embed);
        G.coords.bl_fourvel_to_native(Xnative, ucon_embed, ucon_native);
        DLOOP1 ucov_embed[mu] = 0.;
        DLOOP2 ucov_embed[mu] += gcov_bl[mu][nu] * ucon_embed[nu]; // lower

        // make magnetic fields into a four-vector
        Real B_embed[NVEC], B_native[NVEC], bcon_embed[GR_DIM], bcon_native[GR_DIM];
        B_embed[0] = B(itemp, 0);
        B_embed[1] = B(itemp, 1) / r; // TODO: what to do here? should I reconstruct all vel and B components?
        B_embed[2] = B(itemp, 2) / (r * m::sin(th));
        bcon_embed[0] = 0;
        VLOOP bcon_embed[0] += B_embed[v] * ucov_embed[v+1];
        VLOOP bcon_embed[v+1] = (B_embed[v] + bcon_embed[0] * ucon_embed[v+1]) / ucon_embed[0];

        // convert into native coordinates
        coords.con_vec_to_native(Xnative, bcon_embed, bcon_native);
        Real bcov_native[GR_DIM], bcov_embed[GR_DIM];
        DLOOP1 bcov_embed[mu] = 0.;
        DLOOP2 bcov_embed[mu] += gcov_bl[mu][nu] * bcon_embed[nu]; // lower
        G.lower(bcon_native, bcov_native, k, j, i, Loci::center);

        // convert the 4-vector into 3-vector
        VLOOP B_native[v] = bcon_native[v+1] * ucon_native[0] - bcon_native[0] * ucon_native[v+1];
        
        // save it into B_Save in conserved quantity
        if (i==341 && j==36 && k == 36) printf("HYERIN: i=%d data = (%.3g, %.3g, %.3g) bcon_embed (%.3g, %.3g, %.3g, %.3g) bsq %.3g and native (%.3g, %.3g, %.3g) bsq %.3g \n",
                                    i, B_embed[0], B_embed[1], B_embed[2],
                                    bcon_embed[0], bcon_embed[1], bcon_embed[2], bcon_embed[3], dot(bcon_embed, bcov_embed),
                                    B_native[0], B_native[1], B_native[2], dot(bcon_native, bcov_native));
        VLOOP B_save(v, k, j, i) = B_native[v];
    }

}
