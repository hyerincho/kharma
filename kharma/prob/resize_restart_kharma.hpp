// Load the grid variables up with primitives from an old KHARMA run
#pragma once

#include "decs.hpp"

#include "mesh/mesh.hpp"

// added by Hyerin (10/07/22)
#include "bondi.hpp"
#include "b_flux_ct.hpp"

void ReadFillFile(int i, std::unique_ptr<ParameterInput>& pin);
/**
 * Read the header of an KHARMA HDF5 restart file, and set appropriate parameters
 * Call this before mesh creation!
 */
void ReadKharmaRestartHeader(std::string fname, std::unique_ptr<ParameterInput>& pin);

/**
 * Read data from an KHARMA restart file. Does not support >1 meshblock in Parthenon
 * 
 * Returns stop time tf of the original simulation, for e.g. replicating regression tests
 */
TaskStatus ReadKharmaRestart(MeshBlockData<Real> *rc, ParameterInput *pin);

// newly added by Hyerin (09/06/22)
TaskStatus SetKharmaRestart(MeshBlockData<Real> *rc, IndexDomain domain=IndexDomain::entire, bool coarse=false);

// Hint form resize.hpp
// TODO: (Hyerin) should I do const for x1, x2, x3, var?
KOKKOS_INLINE_FUNCTION void Xtoindex(const GReal XG[GR_DIM],
                                  //Real *x1, Real *x2, Real *x3,
                                   const GridScalar& x1, const GridScalar& x2, const GridScalar& x3,
                                   const GReal dx[GR_DIM], const hsize_t length[GR_DIM], int& iblock,
                                   int& i, int& j, int& k, GReal del[GR_DIM])
{
    //cout << "Hyerin: entered Xtoindex" <<endl;
    Real dx1, dx2, dx3, dx_sum, dx1_min, dx2_min, dx3_min, dx_sum_min;

    // initialize
    iblock =0;
    i = 0;
    j = 0;
    k = 0;
    dx1_min = (XG[1]-x1(iblock,i))*(XG[1]-x1(iblock,i));
    dx2_min = (XG[2]-x2(iblock,j))*(XG[2]-x2(iblock,j));
    dx3_min = (XG[3]-x3(iblock,k))*(XG[3]-x3(iblock,k));
    dx_sum_min = dx1_min + dx2_min + dx3_min;

    for (int iblocktemp = 0; iblocktemp < length[0]; iblocktemp++) {
        // independently searching for minimum for i,j,k
        for (int itemp = 0; itemp < length[1]; itemp++) {
            dx1 = (XG[1]-x1(iblocktemp,itemp))*(XG[1]-x1(iblocktemp,itemp));
            if (dx1 < dx1_min) {
                dx1_min = dx1;
                i = itemp;
            }
        }
        for (int jtemp = 0; jtemp < length[2]; jtemp++) {
            dx2 = (XG[2]-x2(iblocktemp,jtemp))*(XG[2]-x2(iblocktemp,jtemp));
            if (dx2 < dx2_min) {
                dx2_min = dx2;
                j = jtemp;
            }
        }
        for (int ktemp = 0; ktemp < length[3]; ktemp++) {
            dx3 = (XG[3]-x3(iblocktemp,ktemp))*(XG[3]-x3(iblocktemp,ktemp));
            if (dx3 < dx3_min) {
                dx3_min = dx3;
                k = ktemp;
            }
        }
        dx_sum = (XG[1]-x1(iblocktemp,i))*(XG[1]-x1(iblocktemp,i)) + 
                 (XG[2]-x2(iblocktemp,j))*(XG[2]-x2(iblocktemp,j)) + 
                 (XG[3]-x3(iblocktemp,k))*(XG[3]-x3(iblocktemp,k));
        if (dx_sum < dx_sum_min) {
            dx_sum_min = dx_sum;
            iblock = iblocktemp;
        }
    }

    // now construct del
    del[1] = (XG[1] - x1(iblock,i)) / dx[1];
    del[2] = (XG[2] - x2(iblock,j)) / dx[2];
    del[3] = (XG[3] - x3(iblock,k)) / dx[3];
    if (m::abs(dx1_min/m::pow(XG[1],2.))>1.e-8) printf("Xtoindex: dx2 pretty large = %g at r= %g \n",dx1_min, XG[1]);
    if (m::abs(dx2_min)>1.e-8) printf("Xtoindex: dx2 pretty large = %g at th = %g \n",dx2_min, XG[2]);
    if (m::abs(dx3_min)>1.e-8) printf("Xtoindex: dx2 pretty large = %g at phi = %g \n",dx3_min, XG[3]);
}


// linear interpolation for this case. Idea taken from prob/interpolation.hpp
KOKKOS_INLINE_FUNCTION Real linear_interp_kharma(GReal del[GR_DIM], const hsize_t length[GR_DIM], 
                                                const int iblock, const int i, const int j, const int k, const int v, const GridVector& var)
{
    // WARNING: Hyerin (07/20/23) TODO this only supports for one meshblock splitting. I haven't thought about ghost cells that are physical cells

    Real interp;
    int itemp = i;
    int jtemp = j;
    // For ghost zones, we treat each boundary differently:
    // In X1, only stop at the very last cell, where it does not have i+1 to interpolate with
    if (i > length[1] - 2) { itemp = length[1] - 2; del[1] = 1; }
    // In X2, do the same
    if (j > length[2] - 2) { jtemp = length[2] - 2; del[2] = 1; }
    // k auto-wraps. So do all indices for periodic boxes.

    // interpolate in x1 and x2
    interp = var(v, iblock, k, jtemp    , itemp    )*(1. - del[1])*(1. - del[2]) +
             var(v, iblock, k, jtemp + 1, itemp    )*(1. - del[1])*del[2] +
             var(v, iblock, k, jtemp    , itemp + 1)*del[1]*(1. - del[2]) +
             var(v, iblock, k, jtemp + 1, itemp + 1)*del[1]*del[2];

    // then interpolate in x3 if we need
    if (length[3] > 1) {
        interp = (1. - del[3])*interp +
                del[3]*(var(v, iblock, k + 1, jtemp    , itemp    )*(1. - del[1])*(1. - del[2]) +
                        var(v, iblock, k + 1, jtemp + 1, itemp    )*(1. - del[1])*del[2] +
                        var(v, iblock, k + 1, jtemp    , itemp + 1)*del[1]*(1. - del[2]) +
                        var(v, iblock, k + 1, jtemp + 1, itemp + 1)*del[1]*del[2]);
    }

    return interp;
}

KOKKOS_INLINE_FUNCTION void convert_to_utwiddle(const GRCoordinates& G, const CoordinateEmbedding& coords,
                                           const SphBLCoords& bl,  const SphKSCoords& ks, 
                                           const int& k, const int& j, const int& i, Real ucon_bl[GR_DIM], Real u_prim[NVEC])
{
    GReal Xnative[GR_DIM], Xembed[GR_DIM]; //
    G.coord(k, j, i, Loci::center, Xnative);
    G.coord_embed(k, j, i, Loci::center, Xembed);

    // Set u^t to make u^r a 4-vector
    Real gcov_bl[GR_DIM][GR_DIM];
    bl.gcov_embed(Xembed, gcov_bl);
    set_ut(gcov_bl, ucon_bl);

    // Then transform that 4-vector to KS, then to native
    Real ucon_ks[GR_DIM], ucon_mks[GR_DIM];
    ks.vec_from_bl(Xembed, ucon_bl, ucon_ks);
    coords.con_vec_to_native(Xnative, ucon_ks, ucon_mks);

    Real gcon[GR_DIM][GR_DIM];
    G.gcon(Loci::center, j, i, gcon); //TODO: this causes the memory issue!!
    fourvel_to_prim(gcon, ucon_mks, u_prim);

}

KOKKOS_INLINE_FUNCTION void get_prim_restart_kharma(const GReal X[GR_DIM], const VariablePack<Real>& P, const VarMap& m_p, const GReal dx[GR_DIM], const hsize_t length[GR_DIM],
                    const GridScalar& x1, const GridScalar& x2, const GridScalar& x3, const GridScalar& rho, const GridScalar& u, const GridVector& uvec,
                    const int& k, const int& j, const int& i) 
{
    Real rho_temp, u_temp;
    Real u_prim[NVEC]; //, B_prim[NVEC];
    
    GReal del[GR_DIM]; // not really needed now since I am doing nearest neighbor interpolation
    int iblocktemp, itemp, jtemp, ktemp;
    // Interpolate the value at this location from the global grid
    Xtoindex(X, x1, x2, x3, dx, length, iblocktemp, itemp, jtemp, ktemp, del);
    rho_temp = rho(iblocktemp,ktemp,jtemp,itemp);
    u_temp = u(iblocktemp,ktemp,jtemp,itemp);
    VLOOP u_prim[v] = uvec(v,iblocktemp,ktemp,jtemp,itemp);

    P(m_p.RHO, k, j, i) = rho_temp;
    P(m_p.UU, k, j, i) = u_temp;
    P(m_p.U1, k, j, i) = u_prim[0]; 
    P(m_p.U2, k, j, i) = u_prim[1];
    P(m_p.U3, k, j, i) = u_prim[2];

}

KOKKOS_INLINE_FUNCTION void get_B_restart_kharma(const GReal X[GR_DIM], const VariablePack<Real>& P, const VarMap& m_p, const GReal dx[GR_DIM], const hsize_t length[GR_DIM],
                    const GridScalar& x1, const GridScalar& x2, const GridScalar& x3, const GridVector& B, const GridVector& B_save,
                    const int& k, const int& j, const int& i) 
{
    //Real B_prim[NVEC];
    Real B_cons[NVEC];
    
    GReal del[GR_DIM]; // not really needed now since I am doing nearest neighbor interpolation
    int iblocktemp, itemp, jtemp, ktemp;
    // Interpolate the value at this location from the global grid
    Xtoindex(X, x1, x2, x3, dx, length, iblocktemp, itemp, jtemp, ktemp, del);
    //VLOOP B_cons[v] = B(v,iblocktemp,ktemp,jtemp,itemp);
    VLOOP B_cons[v] = linear_interp_kharma(del, length, iblocktemp, itemp, jtemp, ktemp, v, B);

    B_save(0, k, j, i) = B_cons[0];
    B_save(1, k, j, i) = B_cons[1];
    B_save(2, k, j, i) = B_cons[2];

}
