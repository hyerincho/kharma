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

    const Real rs = pin->GetOrAddReal("bondi", "rs", 8.0);
    const Real ur_frac = pin->GetOrAddReal("bondi", "ur_frac", 1.);
    const Real uphi = pin->GetOrAddReal("bondi", "uphi", 0.); 
    auto b_field_type = pin->GetOrAddString("b_field", "type", "none");

    // Set the innermost radius to apply the initialization
    const Real a = pin->GetReal("coordinates", "a");
    const Real rin_default = 1 + m::sqrt(1 - a*a) + 0.1;
    const Real rin_init = pin->GetOrAddReal("gizmo", "r_in", rin_default);

    auto datfn = pin->GetOrAddString("gizmo", "datfn", "none");

    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("gizmo_dat")))
        pmb->packages.Get("GRMHD")->AddParam<std::string>("gizmo_dat", datfn);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("rin_init")))
        pmb->packages.Get("GRMHD")->AddParam<Real>("rin_init", rin_init);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("rs")))
        pmb->packages.Get("GRMHD")->AddParam<Real>("rs", rs);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("ur_frac")))
        pmb->packages.Get("GRMHD")->AddParam<Real>("ur_frac", ur_frac);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("uphi")))
        pmb->packages.Get("GRMHD")->AddParam<Real>("uphi", uphi);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("b_field_type")))
        pmb->packages.Get("GRMHD")->AddParam<std::string>("b_field_type", b_field_type);

    SetGIZMO(rc, IndexDomain::entire);

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

    const Real rs = pmb->packages.Get("GRMHD")->Param<Real>("rs");
    const Real ur_frac = pmb->packages.Get("GRMHD")->Param<Real>("ur_frac");
    const Real uphi = pmb->packages.Get("GRMHD")->Param<Real>("uphi");
    const Real gam = pmb->packages.Get("GRMHD")->Param<Real>("gamma");
    auto datfn = pmb->packages.Get("GRMHD")->Param<std::string>("gizmo_dat");
    auto rin_init = pmb->packages.Get("GRMHD")->Param<Real>("rin_init");
    auto b_field_type = pmb->packages.Get("GRMHD")->Param<std::string>("b_field_type");
    const bool include_B = (b_field_type != "none");
    GridVector B_Save;
    if (include_B) B_Save = rc->Get("B_Save").data;

    // Just the X1 right boundary
    GRCoordinates G = pmb->coords;

    // Set the Bondi conditions wherever we're asked
    auto bounds = coarse ? pmb->c_cellbounds : pmb->cellbounds;

    const IndexRange ib = bounds.GetBoundsI(domain);
    const IndexRange jb = bounds.GetBoundsJ(domain);
    const IndexRange kb = bounds.GetBoundsK(domain);
    
    // GIZMO shell
    // Read the gizmo data file
    std::string fnstr(datfn.c_str());
    std::string dat_type=fnstr.substr(fnstr.find_last_of(".") + 1);
    const bool use_3d = (dat_type != "txt");

    if (! use_3d) {
        std::cout << "GIZMO dat file is txt" << std::endl;
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

                get_prim_gizmo_shell(G, P, m_p, gam, rin_init, rs, ur_frac, uphi, 
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
        length = hdf5_dim[0];
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

        pmb->par_for("gizmo_shell", kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
            KOKKOS_LAMBDA (const int &k, const int &j, const int &i) {

                get_prim_gizmo_shell_3d(G, P, m_p, gam, rin_init, rs, ur_frac, uphi, 
                                      coord_device, rho_device, T_device, v_device, length, k, j, i);

                if (include_B)
                    get_B_gizmo_shell_3d(G, P, rin_init, coord_device, v_device, B_device, B_Save, length, k, j, i);
            }
        );
    }

    return TaskStatus::complete;
}
