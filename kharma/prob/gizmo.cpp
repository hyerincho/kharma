/* 
 *  File: bondi.cpp
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

#include "gizmo.hpp"

#include "floors.hpp"
#include "flux_functions.hpp"

/**
 * Initialization of domain from output of cosmological simulation code GIZMO
 * Note this requires 
 */
TaskStatus InitializeGIZMO(std::shared_ptr<MeshBlockData<Real>>& rc, ParameterInput *pin)
{
    auto pmb = rc->GetBlockPointer();

    const Real mdot = pin->GetOrAddReal("bondi", "mdot", 1.0);
    const Real rs = pin->GetOrAddReal("bondi", "rs", 8.0);
    const Real r_shell = pin->GetOrAddReal("bondi", "r_shell", 0.); 
    const Real ur_frac = pin->GetOrAddReal("bondi", "ur_frac", 1.); 
    const Real uphi = pin->GetOrAddReal("bondi", "uphi", 0.); 
    const Real vacuum_logrho = pin->GetOrAddReal("bondi", "vacuum_logrho", 0.); 
    const Real vacuum_log_u_over_rho = pin->GetOrAddReal("bondi", "vacuum_log_u_over_rho", 0.); 
    auto b_field_type = pin->GetOrAddString("b_field", "type", "none");
    const bool include_B = (b_field_type != "none");

    // Set the innermost radius to apply the initialization
    const Real a = pin->GetReal("coordinates", "a");
    const Real rin_default = 1 + m::sqrt(1 - a*a) + 0.1;
    const Real rin_init = pin->GetOrAddReal("gizmo", "r_in", rin_default);

    auto datfn = pin->GetOrAddString("gizmo", "datfn", "none");

    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("mdot")))
        pmb->packages.Get("GRMHD")->AddParam<Real>("mdot", mdot);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("rs")))
        pmb->packages.Get("GRMHD")->AddParam<Real>("rs", rs);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("gizmo_dat")))
        pmb->packages.Get("GRMHD")->AddParam<std::string>("gizmo_dat", datfn);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("rin_init")))
        pmb->packages.Get("GRMHD")->AddParam<Real>("rin_init", rin_init);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("r_shell")))
        pmb->packages.Get("GRMHD")->AddParam<Real>("r_shell", r_shell);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("ur_frac")))
        pmb->packages.Get("GRMHD")->AddParam<Real>("ur_frac", ur_frac);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("uphi")))
        pmb->packages.Get("GRMHD")->AddParam<Real>("uphi", uphi);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("vacuum_logrho")))
        pmb->packages.Get("GRMHD")->AddParam<Real>("vacuum_logrho", vacuum_logrho);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("vacuum_log_u_over_rho")))
        pmb->packages.Get("GRMHD")->AddParam<Real>("vacuum_log_u_over_rho", vacuum_log_u_over_rho);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("include_B")))
        pmb->packages.Get("GRMHD")->AddParam<Real>("include_B", include_B);

    // Set the interior domain to the analytic solution to begin
    // This tests that PostInitialize will correctly fill ghost zones with the boundary we set
    SetGIZMO(rc, IndexDomain::interior);

    return TaskStatus::complete;
}

TaskStatus SetGIZMO(std::shared_ptr<MeshBlockData<Real>>& rc, IndexDomain domain, bool coarse)
{
    auto pmb = rc->GetBlockPointer();

    //std::cerr << "GIZMO on domain: " << BoundaryName(domain) << std::endl;
    // Don't apply GIZMO initialization to X1 boundaries
    if (domain == IndexDomain::outer_x1 || domain == IndexDomain::inner_x1) {
        return TaskStatus::complete;
    }

    PackIndexMap prims_map, cons_map;
    auto P = GRMHD::PackMHDPrims(rc.get(), prims_map);
    auto U = GRMHD::PackMHDCons(rc.get(), cons_map);
    const VarMap m_u(cons_map, true), m_p(prims_map, false);

    const Real mdot = pmb->packages.Get("GRMHD")->Param<Real>("mdot");
    const Real rs = pmb->packages.Get("GRMHD")->Param<Real>("rs");
    const Real gam = pmb->packages.Get("GRMHD")->Param<Real>("gamma");

    const EMHD::EMHD_parameters& emhd_params = EMHD::GetEMHDParameters(pmb->packages);

    auto datfn = pmb->packages.Get("GRMHD")->Param<std::string>("gizmo_dat");
    auto rin_init = pmb->packages.Get("GRMHD")->Param<Real>("rin_init");
    const Real r_shell = pmb->packages.Get("GRMHD")->Param<Real>("r_shell");
    const Real ur_frac = pmb->packages.Get("GRMHD")->Param<Real>("ur_frac");
    const Real uphi = pmb->packages.Get("GRMHD")->Param<Real>("uphi");
    Real vacuum_logrho = pmb->packages.Get("GRMHD")->Param<Real>("vacuum_logrho");
    Real vacuum_log_u_over_rho = pmb->packages.Get("GRMHD")->Param<Real>("vacuum_log_u_over_rho");
    const bool include_B = pmb->packages.Get("GRMHD")->Param<Real>("include_B");
    GridVector B_Save;
    if (include_B) B_Save = rc->Get("B_Save").data;

    // Just the X1 right boundary
    GRCoordinates G = pmb->coords;
    CoordinateEmbedding cs = G.coords;

    // Set the Bondi conditions wherever we're asked
    auto bounds = coarse ? pmb->c_cellbounds : pmb->cellbounds;

    const IndexRange ib = bounds.GetBoundsI(domain);
    const IndexRange jb = bounds.GetBoundsJ(domain);
    const IndexRange kb = bounds.GetBoundsK(domain);

    std::string fnstr(datfn.c_str());
    std::string dat_type=fnstr.substr(fnstr.find_last_of(".") + 1);
    const bool use_3d = (dat_type != "txt");
    
    // txt file is 1D data
    if (! use_3d) { 
        // Read the gizmo data file
        FILE *fptr = fopen(datfn.c_str(),"r");
        const int datlen = 100000;
        Real *rarr = new double[datlen];
        Real *rhoarr = new double[datlen]; 
        Real *Tarr = new double[datlen]; 
        Real *vrarr = new double[datlen]; 
        Real *Mencarr = new double[datlen]; 
        int length=0, itemp=0;
        while (fscanf(fptr,"%lf %lf %lf %lf %lf\n", &(rarr[itemp]), &(rhoarr[itemp]), &(Tarr[itemp]), &(vrarr[itemp]), &(Mencarr[itemp])) == 5) { // assign the read value to variable, and enter it in array
                itemp++;
        }
        fclose(fptr);
        length = itemp;

        GridVector r_device("r_device", length); 
        GridVector rho_device("rho_device", length); 
        GridVector T_device("T_device", length); 
        GridVector vr_device("vr_device", length); 
        auto r_host = r_device.GetHostMirror();
        auto rho_host = rho_device.GetHostMirror();
        auto T_host = T_device.GetHostMirror();
        auto vr_host = vr_device.GetHostMirror();
        for (itemp = 0; itemp < length; itemp++) {
            r_host(itemp) = rarr[itemp];
            rho_host(itemp) = rhoarr[itemp];
            T_host(itemp) = Tarr[itemp];
            vr_host(itemp) = vrarr[itemp];
        }
        r_device.DeepCopy(r_host);
        rho_device.DeepCopy(rho_host);
        T_device.DeepCopy(T_host);
        vr_device.DeepCopy(vr_host);
            
        Kokkos::fence();

        pmb->par_for("gizmo_shell", kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
            KOKKOS_LAMBDA (const int &k, const int &j, const int &i) {
                // same vacuum conditions at rin_init
                GReal Xshell[GR_DIM] = {0, rin_init, 0, 0};
                int i_sh;
                GReal del_sh;
                XtoindexGIZMO(Xshell, r_device, length, i_sh, del_sh);
                Real vacuum_rho, vacuum_u_over_rho;
                vacuum_rho = rho_device(i_sh)*(1.-del_sh)+rho_device(i_sh+1)*del_sh;
                vacuum_u_over_rho = (T_device(i_sh)*(1.-del_sh)+T_device(i_sh+1)*del_sh)/(gam-1.);

                get_prim_gizmo_shell(G, cs, P, m_p, gam, rin_init, rs, vacuum_rho, vacuum_u_over_rho, 
                    r_device, rho_device, T_device, vr_device, length, k, j, i);
            }
        );
    } else { 
        // hdf5 files are 3D data
        std::cout << "GIZMO dat file is hdf5" << std::endl;

        // read from hdf5 file
        hdf5_open(datfn.c_str());
        hdf5_set_directory("/");
        hsize_t length;
        hsize_t* hdf5_dim;
        hdf5_dim = hdf5_get_attribute_info("PartType0_dimless/Coordinates");
        //printf("HYERIN: the dimensions are (%d, %d)\n", hdf5_dim[0], hdf5_dim[1]);
        length = hdf5_dim[0]; //7372800;
        printf("HYERIN: GIZMO data length is %d \n", length);
        Real *coordarr = new double[length*3];
        Real *rhoarr = new double[length];
        Real *Tarr = new double[length];
        Real *varr = new double[length*3];
        Real *Barr = new double[length*3];
        static hsize_t fdims_vec[] = {length, 3};
        static hsize_t fdims_scl[] = {length};
        hsize_t fstart_vec[] = {0, 0};
        hsize_t fstart_scl[] = {0};
        hdf5_read_array(coordarr, "PartType0_dimless/Coordinates", 2, fdims_vec, fstart_vec,fdims_vec,fdims_vec,fstart_vec,H5T_IEEE_F64LE);
        hdf5_read_array(rhoarr, "PartType0_dimless/Density", 1, fdims_scl, fstart_scl,fdims_scl,fdims_scl,fstart_scl,H5T_IEEE_F64LE);
        hdf5_read_array(Tarr, "PartType0_dimless/Temperature", 1, fdims_scl, fstart_scl,fdims_scl,fdims_scl,fstart_scl,H5T_IEEE_F64LE);
        hdf5_read_array(varr, "PartType0_dimless/Velocities", 2, fdims_vec, fstart_vec,fdims_vec,fdims_vec,fstart_vec,H5T_IEEE_F64LE);
        if (include_B) hdf5_read_array(Barr, "PartType0_dimless/MagneticFields", 2, fdims_vec, fstart_vec,fdims_vec,fdims_vec,fstart_vec,H5T_IEEE_F64LE);
        hdf5_close();
        

        // save in a device arrays
        GridVector coord_device("coord_device", length, 3); 
        GridScalar rho_device("rho_device", length); 
        GridScalar T_device("rho_device", length); 
        GridVector v_device("v_device", length, 3); 
        GridVector B_device("B_device", length, 3); 
        auto coord_host = coord_device.GetHostMirror();
        auto rho_host = rho_device.GetHostMirror();
        auto T_host = T_device.GetHostMirror();
        auto v_host = v_device.GetHostMirror();
        auto B_host = B_device.GetHostMirror();
        int vector_file_index;
        for (int itemp = 0; itemp < length; itemp++) {
            for (int ltemp = 0; ltemp < 3; ltemp++) {
                vector_file_index = 3*itemp+ltemp;
                coord_host(itemp,ltemp) = coordarr[vector_file_index];
                v_host(itemp,ltemp) = varr[vector_file_index];
                if (include_B) B_host(itemp,ltemp) = Barr[vector_file_index];
            }
            rho_host(itemp) = rhoarr[itemp];
            T_host(itemp) = Tarr[itemp];
        }
        coord_device.DeepCopy(coord_host);
        rho_device.DeepCopy(rho_host);
        T_device.DeepCopy(T_host);
        v_device.DeepCopy(v_host);
        if (include_B) B_device.DeepCopy(B_host);
            
        Kokkos::fence();

        int i_sh=0; // TODO! Ask Ben

        pmb->par_for("gizmo_shell", kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
            KOKKOS_LAMBDA (const int &k, const int &j, const int &i) {
                get_prim_gizmo_shell_3d(G, cs, P, m_p, gam, mdot, rs, r_shell, ur_frac, uphi, vacuum_logrho, vacuum_log_u_over_rho, 
                                      coord_device, rho_device, T_device, v_device, length, k, j, i);
                // TODO all flux
                GRMHD::p_to_u(G, P, m_p, gam, k, j, i, U, m_u);

                if (include_B)
                    get_B_gizmo_shell_3d(G, cs, P, r_shell, coord_device, v_device, B_device, B_Save, length, k, j, i);
            }
        );
    }

    return TaskStatus::complete;
}
