/* 
 *  File: resize_restart_kharma.cpp
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

#include "resize_restart_kharma.hpp"

#include "boundaries.hpp"
#include "hdf5_utils.h"
#include "types.hpp"

#include <sys/stat.h>
#include <ctype.h>

// Reads in KHARMA restart file but at a different simulation size

void ReadFillFile(int i, ParameterInput *pin) {
    auto fname_fill = pin->GetOrAddString("resize_restart", "fname_fill" + std::to_string(i), "none");

    if (!(fname_fill == "none")) {
        std::unique_ptr<RestartReader> restartReader;
        restartReader = std::make_unique<RestartReader>(fname_fill.c_str());

        // Load input stream
        std::unique_ptr<ParameterInput> fpinput;
        fpinput = std::make_unique<ParameterInput>();
        auto inputString = restartReader->GetAttr<std::string>("Input", "File");
        std::istringstream is(inputString);
        fpinput->LoadFromStream(is);

        Real fx1min = fpinput->GetReal("parthenon/mesh", "x1min");
        Real fx1max = fpinput->GetReal("parthenon/mesh", "x1max");
        Real fnx1 = fpinput->GetInteger("parthenon/mesh", "nx1");
        Real fnx2 = fpinput->GetInteger("parthenon/mesh", "nx2");
        Real fnx3 = fpinput->GetInteger("parthenon/mesh", "nx3");
        Real fmbnx1 = fpinput->GetInteger("parthenon/meshblock", "nx1");
        Real fmbnx2 = fpinput->GetInteger("parthenon/meshblock", "nx2");
        Real fmbnx3 = fpinput->GetInteger("parthenon/meshblock", "nx3");
        
        restartReader = nullptr;

        pin->SetReal("parthenon/mesh", "restart_x1min_f" + std::to_string(i), fx1min);
        pin->SetReal("parthenon/mesh", "restart_x1max_f" + std::to_string(i), fx1max);
        pin->SetInteger("parthenon/mesh", "restart" + std::to_string(i) + "_nx1", fnx1);
        pin->SetInteger("parthenon/meshblock", "restart" + std::to_string(i) + "_nx1", fmbnx1);
        pin->SetInteger("parthenon/mesh", "restart" + std::to_string(i) + "_nx2", fnx2);
        pin->SetInteger("parthenon/meshblock", "restart" + std::to_string(i) + "_nx2", fmbnx2);
        pin->SetInteger("parthenon/mesh", "restart" + std::to_string(i) + "_nx3", fnx3);
        pin->SetInteger("parthenon/meshblock", "restart" + std::to_string(i) + "_nx3", fmbnx3);
    }
}

void ReadKharmaRestartHeader(std::string fname, ParameterInput *pin)
{
    bool use_dt = pin->GetOrAddBoolean("resize_restart", "use_dt", true);
    bool use_tf = pin->GetOrAddBoolean("resize_restart", "use_tf", false);

    // Read input from restart file 
    // (from external/parthenon/src/parthenon_manager.cpp)
    auto restartReader = std::make_unique<RestartReader>(fname.c_str());

    // Load input stream
    std::unique_ptr<ParameterInput> fpinput;
    fpinput = std::make_unique<ParameterInput>();
    auto inputString = restartReader->GetAttr<std::string>("Input", "File");
    std::istringstream is(inputString);
    fpinput->LoadFromStream(is);

    // TODO(BSP) is there a way to copy all parameters finput->pin and fine-tune later?
    int fnx1, fnx2, fnx3, fmbnx1, fmbnx2, fmbnx3;
    fnx1 = fpinput->GetInteger("parthenon/mesh", "nx1");
    fnx2 = fpinput->GetInteger("parthenon/mesh", "nx2");
    fnx3 = fpinput->GetInteger("parthenon/mesh", "nx3");
    fmbnx1 = fpinput->GetInteger("parthenon/meshblock", "nx1");
    fmbnx2 = fpinput->GetInteger("parthenon/meshblock", "nx2");
    fmbnx3 = fpinput->GetInteger("parthenon/meshblock", "nx3");
    Real fx1min = fpinput->GetReal("parthenon/mesh", "x1min");
    Real fx1max = fpinput->GetReal("parthenon/mesh", "x1max");
    bool fghostzones = fpinput->GetBoolean("parthenon/output1", "ghost_zones");
    int fnghost = fpinput->GetInteger("parthenon/mesh", "nghost");
    auto fBfield = fpinput->GetOrAddString("b_field", "type", "none");
    if (pin->GetOrAddBoolean("resize_restart", "use_restart_size", false)) {
        // This locks the mesh size to be zone-for-zone the same as the iharm3d dump file
        pin->SetInteger("parthenon/mesh", "nx1", fnx1);
        pin->SetInteger("parthenon/mesh", "nx2", fnx2);
        pin->SetInteger("parthenon/mesh", "nx3", fnx3);
        pin->SetInteger("parthenon/meshblock", "nx1", fmbnx1);
        pin->SetInteger("parthenon/meshblock", "nx2", fmbnx2);
        pin->SetInteger("parthenon/meshblock", "nx3", fmbnx3);
    }
    // Record the old values in any case
    pin->SetInteger("parthenon/mesh", "restart_nx1", fnx1);
    pin->SetInteger("parthenon/mesh", "restart_nx2", fnx2);
    pin->SetInteger("parthenon/mesh", "restart_nx3", fnx3);
    pin->SetInteger("parthenon/meshblock", "restart_nx1", fmbnx1);
    pin->SetInteger("parthenon/meshblock", "restart_nx2", fmbnx2);
    pin->SetInteger("parthenon/meshblock", "restart_nx3", fmbnx3);
    pin->SetReal("parthenon/mesh", "restart_x1min", fx1min);
    pin->SetReal("parthenon/mesh", "restart_x1max", fx1max);
    pin->SetInteger("parthenon/mesh", "restart_nghost", fnghost);
    pin->SetBoolean("parthenon/mesh", "restart_ghostzones", fghostzones);
    pin->SetString("b_field", "type", fBfield); // (12/07/22) Hyerin need to test

    Real tNow, dt, tf;
    tNow = restartReader->GetAttr<Real>("Info", "Time");
    dt = restartReader->GetAttr<Real>("Info", "dt");
    tf = fpinput->GetReal("parthenon/time", "tlim");
    //int ncycle = restartReader->GetAttr<int>("Info", "NCycle");

    pin->SetReal("parthenon/time", "start_time", tNow);
    if (use_dt) {
        pin->SetReal("parthenon/time", "dt", dt);
    }
    if (use_tf) {
        pin->SetReal("parthenon/time", "tlim", tf);
    }
    //pin->SetInteger("parthenon/time", "ncycle", ncycle);
    // TODO NSTEP, next tdump/tlog, etc?

    GReal a = fpinput->GetReal("coordinates", "a");
    pin->SetReal("coordinates", "a", a);
    if (fpinput->DoesParameterExist("coordinates", "hslope")) {
        GReal hslope = fpinput->GetReal("coordinates", "hslope");
        pin->SetReal("coordinates", "hslope", hslope);
    }

    // File closed here when restartReader falls out of scope
    restartReader = nullptr;

    for (int i_f = 1; i_f < 7; i_f++) ReadFillFile(i_f, pin);
}

TaskStatus ReadKharmaRestart(std::shared_ptr<MeshBlockData<Real>> rc, ParameterInput *pin)
{
    auto pmb = rc->GetBlockPointer();
    const int num_files_max = 8;

    const int n1tot = pin->GetInteger("parthenon/mesh", "restart_nx1");
    const int n2tot = pin->GetInteger("parthenon/mesh", "restart_nx2");
    const int n3tot = pin->GetInteger("parthenon/mesh", "restart_nx3");
    const hsize_t n1mb = pin->GetInteger("parthenon/meshblock", "restart_nx1");
    const hsize_t n2mb = pin->GetInteger("parthenon/meshblock", "restart_nx2");
    const hsize_t n3mb = pin->GetInteger("parthenon/meshblock", "restart_nx3");
    // file names
    std::vector<std::string> fname_arr;
    for (int i_f = 0; i_f < num_files_max; i_f++) {
        fname_arr.push_back(RestartFileName(i_f, pin));
    }
    //auto fname = pin->GetString("resize_restart", "fname"); // Require this, don't guess. 1st priority
    //auto fname_fill1 = pin->GetOrAddString("resize_restart", "fname_fill1", "none"); // 2nd
    //auto fname_fill2 = pin->GetOrAddString("resize_restart", "fname_fill2", "none"); // 3rd
    //auto fname_fill3 = pin->GetOrAddString("resize_restart", "fname_fill3", "none"); // 4th
    //auto fname_fill4 = pin->GetOrAddString("resize_restart", "fname_fill4", "none"); // 5th
    //auto fname_fill5 = pin->GetOrAddString("resize_restart", "fname_fill5", "none"); // 6th
    //auto fname_fill6 = pin->GetOrAddString("resize_restart", "fname_fill6", "none"); // 7th
    //auto fname_fill7 = pin->GetOrAddString("resize_restart", "fname_fill7", "none"); // 8th
    bool should_fill_arr[num_files_max];
    for (int i_f = 0; i_f < num_files_max; i_f++) { // fill files
        should_fill_arr[i_f] = (!(fname_arr[i_f] == "none"));
    }
    // restart files dimensions
    // TODO Hyerin (07/23/23) assume all the "fill" files have the same dimension.
    const int f_n1tot = pin->GetOrAddInteger("parthenon/mesh", "restart1_nx1", -1);
    const int f_n2tot = pin->GetOrAddInteger("parthenon/mesh", "restart1_nx2", -1);
    const int f_n3tot = pin->GetOrAddInteger("parthenon/mesh", "restart1_nx3", -1);
    const hsize_t f_n1mb = pin->GetOrAddInteger("parthenon/meshblock", "restart1_nx1", -1);
    const hsize_t f_n2mb = pin->GetOrAddInteger("parthenon/meshblock", "restart1_nx2", -1);
    const hsize_t f_n3mb = pin->GetOrAddInteger("parthenon/meshblock", "restart1_nx3", -1);
    // fx1mins, fx1maxs
    Real fx1min_arr[8], fx1max_arr[8];
    fx1min_arr[0] = (pin->GetReal("parthenon/mesh", "restart_x1min"));
    fx1max_arr[0] = (pin->GetReal("parthenon/mesh", "restart_x1max"));
    for (int i_f = 1; i_f < 7; i_f++) { // fill files
        fx1min_arr[i_f] = (pin->GetOrAddReal("parthenon/mesh", "restart_x1min_f" + std::to_string(i_f), -1));
        fx1max_arr[i_f] = (pin->GetOrAddReal("parthenon/mesh", "restart_x1max_f" + std::to_string(i_f), -1));
    }
    //const Real fx1min_f1 = pin->GetOrAddReal("parthenon/mesh", "restart_x1min_f1", -1);
    //const Real fx1max_f1 = pin->GetOrAddReal("parthenon/mesh", "restart_x1max_f1", -1);
    //const Real fx1min_f2 = pin->GetOrAddReal("parthenon/mesh", "restart_x1min_f2", -1);
    //const Real fx1max_f2 = pin->GetOrAddReal("parthenon/mesh", "restart_x1max_f2", -1);
    //const Real fx1min_f3 = pin->GetOrAddReal("parthenon/mesh", "restart_x1min_f3", -1);
    //const Real fx1max_f3 = pin->GetOrAddReal("parthenon/mesh", "restart_x1max_f3", -1);
    //const Real fx1min_f4 = pin->GetOrAddReal("parthenon/mesh", "restart_x1min_f4", -1);
    //const Real fx1max_f4 = pin->GetOrAddReal("parthenon/mesh", "restart_x1max_f4", -1);
    //const Real fx1min_f5 = pin->GetOrAddReal("parthenon/mesh", "restart_x1min_f5", -1);
    //const Real fx1max_f5 = pin->GetOrAddReal("parthenon/mesh", "restart_x1max_f5", -1);
    //const Real fx1min_f6 = pin->GetOrAddReal("parthenon/mesh", "restart_x1min_f6", -1);
    //const Real fx1max_f6 = pin->GetOrAddReal("parthenon/mesh", "restart_x1max_f6", -1);
    const Real x2min = pin->GetReal("parthenon/mesh", "x2min");
    const Real x2max = pin->GetReal("parthenon/mesh", "x2max");
    const Real x3min = pin->GetReal("parthenon/mesh", "x3min");
    const Real x3max = pin->GetReal("parthenon/mesh", "x3max");
    const bool fghostzones = pin->GetBoolean("parthenon/mesh", "restart_ghostzones");
    int fnghost = pin->GetReal("parthenon/mesh", "restart_nghost");
    auto b_field_type = pin->GetOrAddString("b_field", "type", "none");
    int verbose = pin->GetOrAddInteger("debug", "verbose", 0);

    // Derived parameters
    hsize_t nBlocks = (int) (n1tot*n2tot*n3tot)/(n1mb*n2mb*n3mb);
    hsize_t f_nBlocks = (int) (f_n1tot*n2tot*n3tot)/(f_n1mb*n2mb*n3mb);
    const GReal dx[GR_DIM] = {0., (fx1max_arr[0] - fx1min_arr[0])/n1tot,
                                  (x2max - x2min)/n2tot,
                                  (x3max - x3min)/n3tot};
    const bool include_B = (b_field_type != "none");
    auto pkgs = pmb->packages.AllPackages();
    const bool b_ct = (pkgs.count("B_CT"));
    // A placeholder to save the B fields for SeedBField
    GridVector B_Save;
    if (include_B) B_Save = rc->Get("B_Save").data;

    auto& G = pmb->coords;

    // read from a restart file and save it to static GridScalar

    if (!fghostzones) fnghost=0; // reset to 0
    int x3factor=1;
    if (n3tot <= 1) x3factor=0; // if less than 3D, do not add ghosts in x3
    hsize_t length[GR_DIM] = {nBlocks,
                                n1mb + 2 * fnghost,
                                n2mb + 2 * fnghost,
                                n3mb + 2 * fnghost * x3factor}; 
    hsize_t f_length[GR_DIM] = {f_nBlocks,
                                f_n1mb + 2 * fnghost,
                                n2mb   + 2 * fnghost,
                                n3mb   + 2 * fnghost * x3factor}; 
    hsize_t lengthB[GR_DIM], f_lengthB[GR_DIM];
    DLOOP1 lengthB[mu] = length[mu] + (mu > 0? b_ct : 0);
    DLOOP1 f_lengthB[mu] = f_length[mu] + (mu > 0? b_ct : 0);
    const int block_sz = length[0]*length[1]*length[2]*length[3];
    const int blockB_sz = lengthB[0]*lengthB[1]*lengthB[2]*lengthB[3];
    const int f_block_sz = f_length[0]*f_length[1]*f_length[2]*f_length[3];
    const int f_blockB_sz = f_lengthB[0]*f_lengthB[1]*f_lengthB[2]*f_lengthB[3];

    if (MPIRank0() && verbose > 0) {
        std::cout << "Reading mesh size " << n1tot << "x" << n2tot << "x" << n3tot <<
                        " block size " << n1mb << "x" << n2mb << "x" << n3mb << std::endl;
        std::cout << "Reading " << length[0] << " meshblocks of total size " <<
                     length[1] << "x" <<  length[2]<< "x" << length[3] << std::endl;
    }
    
    // read from file and stored in device Hyerin (10/18/2022)
    GridScalar x1_f_device("x1_f_device", length[0], length[1]); 
    GridScalar x2_f_device("x2_f_device", length[0], length[2]); 
    GridScalar x3_f_device("x3_f_device", length[0], length[3]); 
    GridScalar x1B_f_device("x1B_f_device", lengthB[0], lengthB[1]); // For the B fields. If b_ct, they take face values
    GridScalar x2B_f_device("x2B_f_device", lengthB[0], lengthB[2]); 
    GridScalar x3B_f_device("x3B_f_device", lengthB[0], lengthB[3]); 
    GridScalar rho_f_device("rho_f_device", length[0], length[3], length[2], length[1]); 
    GridScalar u_f_device("u_f_device", length[0], length[3], length[2], length[1]); 
    GridVector uvec_f_device("uvec_f_device", NVEC, length[0], length[3], length[2], length[1]); 
    GridVector B_f_device("B_f_device", NVEC, lengthB[0], lengthB[3], lengthB[2], lengthB[1]);
    auto x1_f_host = x1_f_device.GetHostMirror();
    auto x2_f_host = x2_f_device.GetHostMirror();
    auto x3_f_host = x3_f_device.GetHostMirror();
    auto x1B_f_host = x1B_f_device.GetHostMirror();
    auto x2B_f_host = x2B_f_device.GetHostMirror();
    auto x3B_f_host = x3B_f_device.GetHostMirror();
    auto rho_f_host = rho_f_device.GetHostMirror();
    auto u_f_host = u_f_device.GetHostMirror();
    auto uvec_f_host = uvec_f_device.GetHostMirror();
    auto B_f_host = B_f_device.GetHostMirror();
    // Hyerin (09/19/2022) : new attempt to read the file 
    hdf5_open(fname_arr[0].c_str());
    hdf5_set_directory("/");
    Real *rho_file = new double[block_sz];
    Real *u_file = new double[block_sz];
    Real *uvec_file = new double[NVEC*block_sz];
    Real *B_file = new double[NVEC*blockB_sz];
    Real *x1_file = new double[length[0]*length[1]];
    Real *x2_file = new double[length[0]*length[2]];
    Real *x3_file = new double[length[0]*length[3]];
    Real *x1B_file = new double[lengthB[0]*(lengthB[1])];
    Real *x2B_file = new double[lengthB[0]*(lengthB[2])];
    Real *x3B_file = new double[lengthB[0]*(lengthB[3])];
    //static hsize_t fdims[] = {length[0], 1, length[3], length[2], length[1],1}; //outdated
    static hsize_t fdims[] = {length[0], length[3], length[2], length[1]};
    //static hsize_t fdims_vec[] = {length[0], length[3], length[2], length[1],3}; //outdated
    static hsize_t fdims_vec[] = {length[0], NVEC, length[3], length[2], length[1]};
    static hsize_t fdimsB_vec[] = {lengthB[0], NVEC, lengthB[3], lengthB[2], lengthB[1]};
    static hsize_t fdims_x1[] = {length[0], length[1]};
    static hsize_t fdims_x2[] = {length[0], length[2]};
    static hsize_t fdims_x3[] = {length[0], length[3]};
    static hsize_t fdimsB_x1[] = {lengthB[0], lengthB[1]};
    static hsize_t fdimsB_x2[] = {lengthB[0], lengthB[2]};
    static hsize_t fdimsB_x3[] = {lengthB[0], lengthB[3]};
    hsize_t fstart[] = {0, 0, 0, 0};
    hsize_t fstart_vec[] = {0, 0, 0, 0, 0};
    hsize_t fstart_x[] = {0, 0};
    hdf5_read_array(rho_file, "prims.rho", 4, fdims, fstart, fdims, fdims, fstart, H5T_IEEE_F64LE);
    hdf5_read_array(u_file, "prims.u", 4, fdims, fstart, fdims, fdims, fstart, H5T_IEEE_F64LE);
    hdf5_read_array(uvec_file, "prims.uvec", 5, fdims_vec, fstart_vec, fdims_vec, fdims_vec, fstart_vec, H5T_IEEE_F64LE);
    //if (include_B) hdf5_read_array(B_file, "prims.B", 5, fdims_vec, fstart_vec, fdims_vec, fdims_vec, fstart_vec, H5T_IEEE_F64LE);
    if (include_B) {
        if (b_ct) hdf5_read_array(B_file, "cons.fB", 5, fdimsB_vec, fstart_vec, fdimsB_vec, fdimsB_vec, fstart_vec, H5T_IEEE_F64LE);
        else if (pkgs.count("B_FluxCT")) hdf5_read_array(B_file, "cons.B", 5, fdims_vec, fstart_vec, fdims_vec, fdims_vec, fstart_vec, H5T_IEEE_F64LE);
    }
    hdf5_read_array(x1_file, "VolumeLocations/x", 2, fdims_x1, fstart_x, fdims_x1, fdims_x1, fstart_x, H5T_IEEE_F64LE);
    hdf5_read_array(x2_file, "VolumeLocations/y", 2, fdims_x2, fstart_x, fdims_x2, fdims_x2, fstart_x, H5T_IEEE_F64LE);
    hdf5_read_array(x3_file, "VolumeLocations/z", 2, fdims_x3, fstart_x, fdims_x3, fdims_x3, fstart_x, H5T_IEEE_F64LE);
    if (b_ct) {
        hdf5_read_array(x1B_file, "Locations/x", 2, fdimsB_x1, fstart_x, fdimsB_x1, fdimsB_x1, fstart_x, H5T_IEEE_F64LE);
        hdf5_read_array(x2B_file, "Locations/y", 2, fdimsB_x2, fstart_x, fdimsB_x2, fdimsB_x2, fstart_x, H5T_IEEE_F64LE);
        hdf5_read_array(x3B_file, "Locations/z", 2, fdimsB_x3, fstart_x, fdimsB_x3, fdimsB_x3, fstart_x, H5T_IEEE_F64LE);
    }
    hdf5_close();
    
    // fill files dimensions
    static hsize_t fill_fdims[] = {f_length[0], f_length[3], f_length[2], f_length[1]};
    static hsize_t fill_fdims_vec[] = {f_length[0], NVEC, f_length[3], f_length[2], f_length[1]};
    static hsize_t fill_fdimsB_vec[] = {f_lengthB[0], NVEC, f_lengthB[3], f_lengthB[2], f_lengthB[1]};
    static hsize_t fill_fdims_x1[] = {f_length[0], f_length[1]};
    static hsize_t fill_fdims_x2[] = {f_length[0], f_length[2]};
    static hsize_t fill_fdims_x3[] = {f_length[0], f_length[3]};
    static hsize_t fill_fdimsB_x1[] = {f_lengthB[0], f_lengthB[1]};
    static hsize_t fill_fdimsB_x2[] = {f_lengthB[0], f_lengthB[2]};
    static hsize_t fill_fdimsB_x3[] = {f_lengthB[0], f_lengthB[3]};

    // fill # 1
    // TODO: make a function that repeats... idk what's the best way to do it yet though
    GridScalar x1_fill1_device("x1_fill1_device", f_length[0], f_length[1]); 
    GridScalar x2_fill1_device("x2_fill1_device", f_length[0], f_length[2]); 
    GridScalar x3_fill1_device("x3_fill1_device", f_length[0], f_length[3]); 
    GridScalar x1B_fill1_device("x1B_fill1_device", f_lengthB[0], f_lengthB[1]); 
    GridScalar x2B_fill1_device("x2B_fill1_device", f_lengthB[0], f_lengthB[2]); 
    GridScalar x3B_fill1_device("x3B_fill1_device", f_lengthB[0], f_lengthB[3]); 
    GridScalar rho_fill1_device("rho_fill1_device", f_length[0], f_length[3], f_length[2], f_length[1]); 
    GridScalar u_fill1_device("u_fill1_device", f_length[0], f_length[3], f_length[2], f_length[1]); 
    GridVector uvec_fill1_device("uvec_fill1_device", NVEC, f_length[0], f_length[3], f_length[2], f_length[1]); 
    GridVector B_fill1_device("B_fill1_device", NVEC, f_lengthB[0], f_lengthB[3], f_lengthB[2], f_lengthB[1]); 
    auto x1_fill1_host = x1_fill1_device.GetHostMirror();
    auto x2_fill1_host = x2_fill1_device.GetHostMirror();
    auto x3_fill1_host = x3_fill1_device.GetHostMirror();
    auto x1B_fill1_host = x1B_fill1_device.GetHostMirror();
    auto x2B_fill1_host = x2B_fill1_device.GetHostMirror();
    auto x3B_fill1_host = x3B_fill1_device.GetHostMirror();
    auto rho_fill1_host = rho_fill1_device.GetHostMirror();
    auto u_fill1_host = u_fill1_device.GetHostMirror();
    auto uvec_fill1_host = uvec_fill1_device.GetHostMirror();
    auto B_fill1_host = B_fill1_device.GetHostMirror();
    Real *rho_filefill1 = new double[f_block_sz];
    Real *u_filefill1 = new double[f_block_sz];
    Real *uvec_filefill1 = new double[f_block_sz*NVEC];
    Real *B_filefill1 = new double[f_blockB_sz*NVEC];
    Real *x1_filefill1 = new double[f_length[0]*f_length[1]];
    Real *x2_filefill1 = new double[f_length[0]*f_length[2]];
    Real *x3_filefill1 = new double[f_length[0]*f_length[3]];
    Real *x1B_filefill1 = new double[f_lengthB[0]*(f_lengthB[1])];
    Real *x2B_filefill1 = new double[f_lengthB[0]*(f_lengthB[2])];
    Real *x3B_filefill1 = new double[f_lengthB[0]*(f_lengthB[3])];
    if (should_fill_arr[1]) { 
        hdf5_open(fname_arr[1].c_str());
        hdf5_set_directory("/");
        hdf5_read_array(rho_filefill1, "prims.rho", 4, fill_fdims, fstart, fill_fdims, fill_fdims, fstart, H5T_IEEE_F64LE);
        hdf5_read_array(u_filefill1, "prims.u", 4, fill_fdims, fstart, fill_fdims, fill_fdims, fstart, H5T_IEEE_F64LE);
        hdf5_read_array(uvec_filefill1, "prims.uvec", 5, fill_fdims_vec, fstart_vec, fill_fdims_vec, fill_fdims_vec, fstart_vec, H5T_IEEE_F64LE);
        if (include_B) {
            if (b_ct) hdf5_read_array(B_filefill1, "cons.fB", 5, fill_fdimsB_vec, fstart_vec, fill_fdimsB_vec, fill_fdimsB_vec, fstart_vec,H5T_IEEE_F64LE);
            else if (pkgs.count("B_FluxCT")) hdf5_read_array(B_filefill1, "cons.B", 5, fill_fdims_vec, fstart_vec, fill_fdims_vec, fill_fdims_vec, fstart_vec,H5T_IEEE_F64LE);
        }
        hdf5_read_array(x1_filefill1, "VolumeLocations/x", 2, fill_fdims_x1, fstart_x, fill_fdims_x1, fill_fdims_x1, fstart_x, H5T_IEEE_F64LE);
        hdf5_read_array(x2_filefill1, "VolumeLocations/y", 2, fill_fdims_x2, fstart_x, fill_fdims_x2, fill_fdims_x2, fstart_x, H5T_IEEE_F64LE);
        hdf5_read_array(x3_filefill1, "VolumeLocations/z", 2, fill_fdims_x3, fstart_x, fill_fdims_x3, fill_fdims_x3, fstart_x, H5T_IEEE_F64LE);
        if (b_ct) {
            hdf5_read_array(x1B_filefill1, "Locations/x", 2, fill_fdimsB_x1, fstart_x, fill_fdimsB_x1, fill_fdimsB_x1, fstart_x, H5T_IEEE_F64LE);
            hdf5_read_array(x2B_filefill1, "Locations/y", 2, fill_fdimsB_x2, fstart_x, fill_fdimsB_x2, fill_fdimsB_x2, fstart_x, H5T_IEEE_F64LE);
            hdf5_read_array(x3B_filefill1, "Locations/z", 2, fill_fdimsB_x3, fstart_x, fill_fdimsB_x3, fill_fdimsB_x3, fstart_x, H5T_IEEE_F64LE);
        }
        hdf5_close();
    }
    
    // fill # 2
    GridScalar x1_fill2_device("x1_fill2_device", f_length[0], f_length[1]); 
    GridScalar x2_fill2_device("x2_fill2_device", f_length[0], f_length[2]); 
    GridScalar x3_fill2_device("x3_fill2_device", f_length[0], f_length[3]); 
    GridScalar x1B_fill2_device("x1B_fill2_device", f_lengthB[0], f_lengthB[1]); 
    GridScalar x2B_fill2_device("x2B_fill2_device", f_lengthB[0], f_lengthB[2]); 
    GridScalar x3B_fill2_device("x3B_fill2_device", f_lengthB[0], f_lengthB[3]); 
    GridScalar rho_fill2_device("rho_fill2_device", f_length[0], f_length[3], f_length[2], f_length[1]); 
    GridScalar u_fill2_device("u_fill2_device", f_length[0], f_length[3], f_length[2], f_length[1]); 
    GridVector uvec_fill2_device("uvec_fill2_device", NVEC, f_length[0], f_length[3], f_length[2], f_length[1]); 
    GridVector B_fill2_device("B_fill2_device", NVEC, f_lengthB[0], f_lengthB[3], f_lengthB[2], f_lengthB[1]); 
    auto x1_fill2_host = x1_fill2_device.GetHostMirror();
    auto x2_fill2_host = x2_fill2_device.GetHostMirror();
    auto x3_fill2_host = x3_fill2_device.GetHostMirror();
    auto x1B_fill2_host = x1B_fill2_device.GetHostMirror();
    auto x2B_fill2_host = x2B_fill2_device.GetHostMirror();
    auto x3B_fill2_host = x3B_fill2_device.GetHostMirror();
    auto rho_fill2_host = rho_fill2_device.GetHostMirror();
    auto u_fill2_host = u_fill2_device.GetHostMirror();
    auto uvec_fill2_host = uvec_fill2_device.GetHostMirror();
    auto B_fill2_host = B_fill2_device.GetHostMirror();
    Real *rho_filefill2 = new double[f_block_sz];
    Real *u_filefill2 = new double[f_block_sz];
    Real *uvec_filefill2 = new double[f_block_sz*NVEC];
    Real *B_filefill2 = new double[f_blockB_sz*NVEC];
    Real *x1_filefill2 = new double[f_length[0]*f_length[1]];
    Real *x2_filefill2 = new double[f_length[0]*f_length[2]];
    Real *x3_filefill2 = new double[f_length[0]*f_length[3]];
    Real *x1B_filefill2 = new double[f_lengthB[0]*(f_lengthB[1])];
    Real *x2B_filefill2 = new double[f_lengthB[0]*(f_lengthB[2])];
    Real *x3B_filefill2 = new double[f_lengthB[0]*(f_lengthB[3])];
    if (should_fill_arr[2]) { 
        hdf5_open(fname_arr[2].c_str());
        hdf5_set_directory("/");
        hdf5_read_array(rho_filefill2, "prims.rho", 4, fill_fdims, fstart, fill_fdims, fill_fdims, fstart, H5T_IEEE_F64LE);
        hdf5_read_array(u_filefill2, "prims.u", 4, fill_fdims, fstart, fill_fdims, fill_fdims, fstart, H5T_IEEE_F64LE);
        hdf5_read_array(uvec_filefill2, "prims.uvec", 5, fill_fdims_vec, fstart_vec, fill_fdims_vec, fill_fdims_vec, fstart_vec, H5T_IEEE_F64LE);
        if (include_B) {
            if (b_ct) hdf5_read_array(B_filefill2, "cons.fB", 5, fill_fdimsB_vec, fstart_vec, fill_fdimsB_vec, fill_fdimsB_vec, fstart_vec,H5T_IEEE_F64LE);
            else if (pkgs.count("B_FluxCT")) hdf5_read_array(B_filefill2, "cons.B", 5, fill_fdims_vec, fstart_vec, fill_fdims_vec, fill_fdims_vec, fstart_vec,H5T_IEEE_F64LE);
        }
        hdf5_read_array(x1_filefill2, "VolumeLocations/x", 2, fill_fdims_x1, fstart_x, fill_fdims_x1, fill_fdims_x1, fstart_x, H5T_IEEE_F64LE);
        hdf5_read_array(x2_filefill2, "VolumeLocations/y", 2, fill_fdims_x2, fstart_x, fill_fdims_x2, fill_fdims_x2, fstart_x, H5T_IEEE_F64LE);
        hdf5_read_array(x3_filefill2, "VolumeLocations/z", 2, fill_fdims_x3, fstart_x, fill_fdims_x3, fill_fdims_x3, fstart_x, H5T_IEEE_F64LE);
        if (b_ct) {
            hdf5_read_array(x1B_filefill2, "Locations/x", 2, fill_fdimsB_x1, fstart_x, fill_fdimsB_x1, fill_fdimsB_x1, fstart_x, H5T_IEEE_F64LE);
            hdf5_read_array(x2B_filefill2, "Locations/y", 2, fill_fdimsB_x2, fstart_x, fill_fdimsB_x2, fill_fdimsB_x2, fstart_x, H5T_IEEE_F64LE);
            hdf5_read_array(x3B_filefill2, "Locations/z", 2, fill_fdimsB_x3, fstart_x, fill_fdimsB_x3, fill_fdimsB_x3, fstart_x, H5T_IEEE_F64LE);
        }
        hdf5_close();
    }
    
    // fill # 3
    GridScalar x1_fill3_device("x1_fill3_device", f_length[0], f_length[1]); 
    GridScalar x2_fill3_device("x2_fill3_device", f_length[0], f_length[2]); 
    GridScalar x3_fill3_device("x3_fill3_device", f_length[0], f_length[3]); 
    GridScalar x1B_fill3_device("x1B_fill3_device", f_lengthB[0], f_lengthB[1]); 
    GridScalar x2B_fill3_device("x2B_fill3_device", f_lengthB[0], f_lengthB[2]); 
    GridScalar x3B_fill3_device("x3B_fill3_device", f_lengthB[0], f_lengthB[3]); 
    GridScalar rho_fill3_device("rho_fill3_device", f_length[0], f_length[3], f_length[2], f_length[1]); 
    GridScalar u_fill3_device("u_fill3_device", f_length[0], f_length[3], f_length[2], f_length[1]); 
    GridVector uvec_fill3_device("uvec_fill3_device", NVEC, f_length[0], f_length[3], f_length[2], f_length[1]); 
    GridVector B_fill3_device("B_fill3_device", NVEC, f_lengthB[0], f_lengthB[3], f_lengthB[2], f_lengthB[1]); 
    auto x1_fill3_host = x1_fill3_device.GetHostMirror();
    auto x2_fill3_host = x2_fill3_device.GetHostMirror();
    auto x3_fill3_host = x3_fill3_device.GetHostMirror();
    auto x1B_fill3_host = x1B_fill3_device.GetHostMirror();
    auto x2B_fill3_host = x2B_fill3_device.GetHostMirror();
    auto x3B_fill3_host = x3B_fill3_device.GetHostMirror();
    auto rho_fill3_host = rho_fill3_device.GetHostMirror();
    auto u_fill3_host = u_fill3_device.GetHostMirror();
    auto uvec_fill3_host = uvec_fill3_device.GetHostMirror();
    auto B_fill3_host = B_fill3_device.GetHostMirror();
    Real *rho_filefill3 = new double[f_block_sz];
    Real *u_filefill3 = new double[f_block_sz];
    Real *uvec_filefill3 = new double[f_block_sz*NVEC];
    Real *B_filefill3 = new double[f_blockB_sz*NVEC];
    Real *x1_filefill3 = new double[f_length[0]*f_length[1]];
    Real *x2_filefill3 = new double[f_length[0]*f_length[2]];
    Real *x3_filefill3 = new double[f_length[0]*f_length[3]];
    Real *x1B_filefill3 = new double[f_lengthB[0]*(f_lengthB[1])];
    Real *x2B_filefill3 = new double[f_lengthB[0]*(f_lengthB[2])];
    Real *x3B_filefill3 = new double[f_lengthB[0]*(f_lengthB[3])];
    if (should_fill_arr[3]) { 
        hdf5_open(fname_arr[3].c_str());
        hdf5_set_directory("/");
        hdf5_read_array(rho_filefill3, "prims.rho", 4, fill_fdims, fstart, fill_fdims, fill_fdims, fstart, H5T_IEEE_F64LE);
        hdf5_read_array(u_filefill3, "prims.u", 4, fill_fdims, fstart, fill_fdims, fill_fdims, fstart, H5T_IEEE_F64LE);
        hdf5_read_array(uvec_filefill3, "prims.uvec", 5, fill_fdims_vec, fstart_vec, fill_fdims_vec, fill_fdims_vec, fstart_vec, H5T_IEEE_F64LE);
        if (include_B) {
            if (b_ct) hdf5_read_array(B_filefill3, "cons.fB", 5, fill_fdimsB_vec, fstart_vec, fill_fdimsB_vec, fill_fdimsB_vec, fstart_vec,H5T_IEEE_F64LE);
            else if (pkgs.count("B_FluxCT")) hdf5_read_array(B_filefill3, "cons.B", 5, fill_fdims_vec, fstart_vec, fill_fdims_vec, fill_fdims_vec, fstart_vec,H5T_IEEE_F64LE);
        }
        hdf5_read_array(x1_filefill3, "VolumeLocations/x", 2, fill_fdims_x1, fstart_x, fill_fdims_x1, fill_fdims_x1, fstart_x, H5T_IEEE_F64LE);
        hdf5_read_array(x2_filefill3, "VolumeLocations/y", 2, fill_fdims_x2, fstart_x, fill_fdims_x2, fill_fdims_x2, fstart_x, H5T_IEEE_F64LE);
        hdf5_read_array(x3_filefill3, "VolumeLocations/z", 2, fill_fdims_x3, fstart_x, fill_fdims_x3, fill_fdims_x3, fstart_x, H5T_IEEE_F64LE);
        if (b_ct) {
            hdf5_read_array(x1B_filefill3, "Locations/x", 2, fill_fdimsB_x1, fstart_x, fill_fdimsB_x1, fill_fdimsB_x1, fstart_x, H5T_IEEE_F64LE);
            hdf5_read_array(x2B_filefill3, "Locations/y", 2, fill_fdimsB_x2, fstart_x, fill_fdimsB_x2, fill_fdimsB_x2, fstart_x, H5T_IEEE_F64LE);
            hdf5_read_array(x3B_filefill3, "Locations/z", 2, fill_fdimsB_x3, fstart_x, fill_fdimsB_x3, fill_fdimsB_x3, fstart_x, H5T_IEEE_F64LE);
        }
        hdf5_close();
    }
    
    // fill # 4
    GridScalar x1_fill4_device("x1_fill4_device", f_length[0], f_length[1]); 
    GridScalar x2_fill4_device("x2_fill4_device", f_length[0], f_length[2]); 
    GridScalar x3_fill4_device("x3_fill4_device", f_length[0], f_length[3]); 
    GridScalar x1B_fill4_device("x1B_fill4_device", f_lengthB[0], f_lengthB[1]); 
    GridScalar x2B_fill4_device("x2B_fill4_device", f_lengthB[0], f_lengthB[2]); 
    GridScalar x3B_fill4_device("x3B_fill4_device", f_lengthB[0], f_lengthB[3]); 
    GridScalar rho_fill4_device("rho_fill4_device", f_length[0], f_length[3], f_length[2], f_length[1]); 
    GridScalar u_fill4_device("u_fill4_device", f_length[0], f_length[3], f_length[2], f_length[1]); 
    GridVector uvec_fill4_device("uvec_fill4_device", NVEC, f_length[0], f_length[3], f_length[2], f_length[1]); 
    GridVector B_fill4_device("B_fill4_device", NVEC, f_lengthB[0], f_lengthB[3], f_lengthB[2], f_lengthB[1]); 
    auto x1_fill4_host = x1_fill4_device.GetHostMirror();
    auto x2_fill4_host = x2_fill4_device.GetHostMirror();
    auto x3_fill4_host = x3_fill4_device.GetHostMirror();
    auto x1B_fill4_host = x1B_fill4_device.GetHostMirror();
    auto x2B_fill4_host = x2B_fill4_device.GetHostMirror();
    auto x3B_fill4_host = x3B_fill4_device.GetHostMirror();
    auto rho_fill4_host = rho_fill4_device.GetHostMirror();
    auto u_fill4_host = u_fill4_device.GetHostMirror();
    auto uvec_fill4_host = uvec_fill4_device.GetHostMirror();
    auto B_fill4_host = B_fill4_device.GetHostMirror();
    Real *rho_filefill4 = new double[f_block_sz];
    Real *u_filefill4 = new double[f_block_sz];
    Real *uvec_filefill4 = new double[f_block_sz*NVEC];
    Real *B_filefill4 = new double[f_blockB_sz*NVEC];
    Real *x1_filefill4 = new double[f_length[0]*f_length[1]];
    Real *x2_filefill4 = new double[f_length[0]*f_length[2]];
    Real *x3_filefill4 = new double[f_length[0]*f_length[3]];
    Real *x1B_filefill4 = new double[f_lengthB[0]*(f_lengthB[1])];
    Real *x2B_filefill4 = new double[f_lengthB[0]*(f_lengthB[2])];
    Real *x3B_filefill4 = new double[f_lengthB[0]*(f_lengthB[3])];
    if (should_fill_arr[4]) { 
        hdf5_open(fname_arr[4].c_str());
        hdf5_set_directory("/");
        hdf5_read_array(rho_filefill4, "prims.rho", 4, fill_fdims, fstart, fill_fdims, fill_fdims, fstart, H5T_IEEE_F64LE);
        hdf5_read_array(u_filefill4, "prims.u", 4, fill_fdims, fstart, fill_fdims, fill_fdims, fstart, H5T_IEEE_F64LE);
        hdf5_read_array(uvec_filefill4, "prims.uvec", 5, fill_fdims_vec, fstart_vec, fill_fdims_vec, fill_fdims_vec, fstart_vec, H5T_IEEE_F64LE);
        if (include_B) {
            if (b_ct) hdf5_read_array(B_filefill4, "cons.fB", 5, fill_fdimsB_vec, fstart_vec, fill_fdimsB_vec, fill_fdimsB_vec, fstart_vec,H5T_IEEE_F64LE);
            else if (pkgs.count("B_FluxCT")) hdf5_read_array(B_filefill4, "cons.B", 5, fill_fdims_vec, fstart_vec, fill_fdims_vec, fill_fdims_vec, fstart_vec,H5T_IEEE_F64LE);
        }
        hdf5_read_array(x1_filefill4, "VolumeLocations/x", 2, fill_fdims_x1, fstart_x, fill_fdims_x1, fill_fdims_x1, fstart_x, H5T_IEEE_F64LE);
        hdf5_read_array(x2_filefill4, "VolumeLocations/y", 2, fill_fdims_x2, fstart_x, fill_fdims_x2, fill_fdims_x2, fstart_x, H5T_IEEE_F64LE);
        hdf5_read_array(x3_filefill4, "VolumeLocations/z", 2, fill_fdims_x3, fstart_x, fill_fdims_x3, fill_fdims_x3, fstart_x, H5T_IEEE_F64LE);
        if (b_ct) {
            hdf5_read_array(x1B_filefill4, "Locations/x", 2, fill_fdimsB_x1, fstart_x, fill_fdimsB_x1, fill_fdimsB_x1, fstart_x, H5T_IEEE_F64LE);
            hdf5_read_array(x2B_filefill4, "Locations/y", 2, fill_fdimsB_x2, fstart_x, fill_fdimsB_x2, fill_fdimsB_x2, fstart_x, H5T_IEEE_F64LE);
            hdf5_read_array(x3B_filefill4, "Locations/z", 2, fill_fdimsB_x3, fstart_x, fill_fdimsB_x3, fill_fdimsB_x3, fstart_x, H5T_IEEE_F64LE);
        }
        hdf5_close();
    }
    
    // fill # 5
    GridScalar x1_fill5_device("x1_fill5_device", f_length[0], f_length[1]); 
    GridScalar x2_fill5_device("x2_fill5_device", f_length[0], f_length[2]); 
    GridScalar x3_fill5_device("x3_fill5_device", f_length[0], f_length[3]); 
    GridScalar x1B_fill5_device("x1B_fill5_device", f_lengthB[0], f_lengthB[1]); 
    GridScalar x2B_fill5_device("x2B_fill5_device", f_lengthB[0], f_lengthB[2]); 
    GridScalar x3B_fill5_device("x3B_fill5_device", f_lengthB[0], f_lengthB[3]); 
    GridScalar rho_fill5_device("rho_fill5_device", f_length[0], f_length[3], f_length[2], f_length[1]); 
    GridScalar u_fill5_device("u_fill5_device", f_length[0], f_length[3], f_length[2], f_length[1]); 
    GridVector uvec_fill5_device("uvec_fill5_device", NVEC, f_length[0], f_length[3], f_length[2], f_length[1]); 
    GridVector B_fill5_device("B_fill5_device", NVEC, f_lengthB[0], f_lengthB[3], f_lengthB[2], f_lengthB[1]); 
    auto x1_fill5_host = x1_fill5_device.GetHostMirror();
    auto x2_fill5_host = x2_fill5_device.GetHostMirror();
    auto x3_fill5_host = x3_fill5_device.GetHostMirror();
    auto x1B_fill5_host = x1B_fill5_device.GetHostMirror();
    auto x2B_fill5_host = x2B_fill5_device.GetHostMirror();
    auto x3B_fill5_host = x3B_fill5_device.GetHostMirror();
    auto rho_fill5_host = rho_fill5_device.GetHostMirror();
    auto u_fill5_host = u_fill5_device.GetHostMirror();
    auto uvec_fill5_host = uvec_fill5_device.GetHostMirror();
    auto B_fill5_host = B_fill5_device.GetHostMirror();
    Real *rho_filefill5 = new double[f_block_sz];
    Real *u_filefill5 = new double[f_block_sz];
    Real *uvec_filefill5 = new double[f_block_sz*NVEC];
    Real *B_filefill5 = new double[f_blockB_sz*NVEC];
    Real *x1_filefill5 = new double[f_length[0]*f_length[1]];
    Real *x2_filefill5 = new double[f_length[0]*f_length[2]];
    Real *x3_filefill5 = new double[f_length[0]*f_length[3]];
    Real *x1B_filefill5 = new double[f_lengthB[0]*(f_lengthB[1])];
    Real *x2B_filefill5 = new double[f_lengthB[0]*(f_lengthB[2])];
    Real *x3B_filefill5 = new double[f_lengthB[0]*(f_lengthB[3])];
    if (should_fill_arr[5]) { 
        hdf5_open(fname_arr[5].c_str());
        hdf5_set_directory("/");
        hdf5_read_array(rho_filefill5, "prims.rho", 4, fill_fdims, fstart, fill_fdims, fill_fdims, fstart, H5T_IEEE_F64LE);
        hdf5_read_array(u_filefill5, "prims.u", 4, fill_fdims, fstart, fill_fdims, fill_fdims, fstart, H5T_IEEE_F64LE);
        hdf5_read_array(uvec_filefill5, "prims.uvec", 5, fill_fdims_vec, fstart_vec, fill_fdims_vec, fill_fdims_vec, fstart_vec, H5T_IEEE_F64LE);
        if (include_B) {
            if (b_ct) hdf5_read_array(B_filefill5, "cons.fB", 5, fill_fdimsB_vec, fstart_vec, fill_fdimsB_vec, fill_fdimsB_vec, fstart_vec,H5T_IEEE_F64LE);
            else if (pkgs.count("B_FluxCT")) hdf5_read_array(B_filefill5, "cons.B", 5, fill_fdims_vec, fstart_vec, fill_fdims_vec, fill_fdims_vec, fstart_vec,H5T_IEEE_F64LE);
        }
        hdf5_read_array(x1_filefill5, "VolumeLocations/x", 2, fill_fdims_x1, fstart_x, fill_fdims_x1, fill_fdims_x1, fstart_x, H5T_IEEE_F64LE);
        hdf5_read_array(x2_filefill5, "VolumeLocations/y", 2, fill_fdims_x2, fstart_x, fill_fdims_x2, fill_fdims_x2, fstart_x, H5T_IEEE_F64LE);
        hdf5_read_array(x3_filefill5, "VolumeLocations/z", 2, fill_fdims_x3, fstart_x, fill_fdims_x3, fill_fdims_x3, fstart_x, H5T_IEEE_F64LE);
        if (b_ct) {
            hdf5_read_array(x1B_filefill5, "Locations/x", 2, fill_fdimsB_x1, fstart_x, fill_fdimsB_x1, fill_fdimsB_x1, fstart_x, H5T_IEEE_F64LE);
            hdf5_read_array(x2B_filefill5, "Locations/y", 2, fill_fdimsB_x2, fstart_x, fill_fdimsB_x2, fill_fdimsB_x2, fstart_x, H5T_IEEE_F64LE);
            hdf5_read_array(x3B_filefill5, "Locations/z", 2, fill_fdimsB_x3, fstart_x, fill_fdimsB_x3, fill_fdimsB_x3, fstart_x, H5T_IEEE_F64LE);
        }
        hdf5_close();
    }
    
    // fill # 6
    GridScalar x1_fill6_device("x1_fill6_device", f_length[0], f_length[1]); 
    GridScalar x2_fill6_device("x2_fill6_device", f_length[0], f_length[2]); 
    GridScalar x3_fill6_device("x3_fill6_device", f_length[0], f_length[3]); 
    GridScalar x1B_fill6_device("x1B_fill6_device", f_lengthB[0], f_lengthB[1]); 
    GridScalar x2B_fill6_device("x2B_fill6_device", f_lengthB[0], f_lengthB[2]); 
    GridScalar x3B_fill6_device("x3B_fill6_device", f_lengthB[0], f_lengthB[3]); 
    GridScalar rho_fill6_device("rho_fill6_device", f_length[0], f_length[3], f_length[2], f_length[1]); 
    GridScalar u_fill6_device("u_fill6_device", f_length[0], f_length[3], f_length[2], f_length[1]); 
    GridVector uvec_fill6_device("uvec_fill6_device", NVEC, f_length[0], f_length[3], f_length[2], f_length[1]); 
    GridVector B_fill6_device("B_fill6_device", NVEC, f_lengthB[0], f_lengthB[3], f_lengthB[2], f_lengthB[1]); 
    auto x1_fill6_host = x1_fill6_device.GetHostMirror();
    auto x2_fill6_host = x2_fill6_device.GetHostMirror();
    auto x3_fill6_host = x3_fill6_device.GetHostMirror();
    auto x1B_fill6_host = x1B_fill6_device.GetHostMirror();
    auto x2B_fill6_host = x2B_fill6_device.GetHostMirror();
    auto x3B_fill6_host = x3B_fill6_device.GetHostMirror();
    auto rho_fill6_host = rho_fill6_device.GetHostMirror();
    auto u_fill6_host = u_fill6_device.GetHostMirror();
    auto uvec_fill6_host = uvec_fill6_device.GetHostMirror();
    auto B_fill6_host = B_fill6_device.GetHostMirror();
    Real *rho_filefill6 = new double[f_block_sz];
    Real *u_filefill6 = new double[f_block_sz];
    Real *uvec_filefill6 = new double[f_block_sz*NVEC];
    Real *B_filefill6 = new double[f_blockB_sz*NVEC];
    Real *x1_filefill6 = new double[f_length[0]*f_length[1]];
    Real *x2_filefill6 = new double[f_length[0]*f_length[2]];
    Real *x3_filefill6 = new double[f_length[0]*f_length[3]];
    Real *x1B_filefill6 = new double[f_lengthB[0]*(f_lengthB[1])];
    Real *x2B_filefill6 = new double[f_lengthB[0]*(f_lengthB[2])];
    Real *x3B_filefill6 = new double[f_lengthB[0]*(f_lengthB[3])];
    if (should_fill_arr[6]) { 
        hdf5_open(fname_arr[6].c_str());
        hdf5_set_directory("/");
        hdf5_read_array(rho_filefill6, "prims.rho", 4, fill_fdims, fstart, fill_fdims, fill_fdims, fstart, H5T_IEEE_F64LE);
        hdf5_read_array(u_filefill6, "prims.u", 4, fill_fdims, fstart, fill_fdims, fill_fdims, fstart, H5T_IEEE_F64LE);
        hdf5_read_array(uvec_filefill6, "prims.uvec", 5, fill_fdims_vec, fstart_vec, fill_fdims_vec, fill_fdims_vec, fstart_vec, H5T_IEEE_F64LE);
        if (include_B) {
            if (b_ct) hdf5_read_array(B_filefill6, "cons.fB", 5, fill_fdimsB_vec, fstart_vec, fill_fdimsB_vec, fill_fdimsB_vec, fstart_vec,H5T_IEEE_F64LE);
            else if (pkgs.count("B_FluxCT")) hdf5_read_array(B_filefill6, "cons.B", 5, fill_fdims_vec, fstart_vec, fill_fdims_vec, fill_fdims_vec, fstart_vec,H5T_IEEE_F64LE);
        }
        hdf5_read_array(x1_filefill6, "VolumeLocations/x", 2, fill_fdims_x1, fstart_x, fill_fdims_x1, fill_fdims_x1, fstart_x, H5T_IEEE_F64LE);
        hdf5_read_array(x2_filefill6, "VolumeLocations/y", 2, fill_fdims_x2, fstart_x, fill_fdims_x2, fill_fdims_x2, fstart_x, H5T_IEEE_F64LE);
        hdf5_read_array(x3_filefill6, "VolumeLocations/z", 2, fill_fdims_x3, fstart_x, fill_fdims_x3, fill_fdims_x3, fstart_x, H5T_IEEE_F64LE);
        if (b_ct) {
            hdf5_read_array(x1B_filefill6, "Locations/x", 2, fill_fdimsB_x1, fstart_x, fill_fdimsB_x1, fill_fdimsB_x1, fstart_x, H5T_IEEE_F64LE);
            hdf5_read_array(x2B_filefill6, "Locations/y", 2, fill_fdimsB_x2, fstart_x, fill_fdimsB_x2, fill_fdimsB_x2, fstart_x, H5T_IEEE_F64LE);
            hdf5_read_array(x3B_filefill6, "Locations/z", 2, fill_fdimsB_x3, fstart_x, fill_fdimsB_x3, fill_fdimsB_x3, fstart_x, H5T_IEEE_F64LE);
        }
        hdf5_close();
    }
    
    // fill # 7
    GridScalar x1_fill7_device("x1_fill7_device", f_length[0], f_length[1]); 
    GridScalar x2_fill7_device("x2_fill7_device", f_length[0], f_length[2]); 
    GridScalar x3_fill7_device("x3_fill7_device", f_length[0], f_length[3]); 
    GridScalar x1B_fill7_device("x1B_fill7_device", f_lengthB[0], f_lengthB[1]); 
    GridScalar x2B_fill7_device("x2B_fill7_device", f_lengthB[0], f_lengthB[2]); 
    GridScalar x3B_fill7_device("x3B_fill7_device", f_lengthB[0], f_lengthB[3]); 
    GridScalar rho_fill7_device("rho_fill7_device", f_length[0], f_length[3], f_length[2], f_length[1]); 
    GridScalar u_fill7_device("u_fill7_device", f_length[0], f_length[3], f_length[2], f_length[1]); 
    GridVector uvec_fill7_device("uvec_fill7_device", NVEC, f_length[0], f_length[3], f_length[2], f_length[1]); 
    GridVector B_fill7_device("B_fill7_device", NVEC, f_lengthB[0], f_lengthB[3], f_lengthB[2], f_lengthB[1]); 
    auto x1_fill7_host = x1_fill7_device.GetHostMirror();
    auto x2_fill7_host = x2_fill7_device.GetHostMirror();
    auto x3_fill7_host = x3_fill7_device.GetHostMirror();
    auto x1B_fill7_host = x1B_fill7_device.GetHostMirror();
    auto x2B_fill7_host = x2B_fill7_device.GetHostMirror();
    auto x3B_fill7_host = x3B_fill7_device.GetHostMirror();
    auto rho_fill7_host = rho_fill7_device.GetHostMirror();
    auto u_fill7_host = u_fill7_device.GetHostMirror();
    auto uvec_fill7_host = uvec_fill7_device.GetHostMirror();
    auto B_fill7_host = B_fill7_device.GetHostMirror();
    Real *rho_filefill7 = new double[f_block_sz];
    Real *u_filefill7 = new double[f_block_sz];
    Real *uvec_filefill7 = new double[f_block_sz*NVEC];
    Real *B_filefill7 = new double[f_blockB_sz*NVEC];
    Real *x1_filefill7 = new double[f_length[0]*f_length[1]];
    Real *x2_filefill7 = new double[f_length[0]*f_length[2]];
    Real *x3_filefill7 = new double[f_length[0]*f_length[3]];
    Real *x1B_filefill7 = new double[f_lengthB[0]*(f_lengthB[1])];
    Real *x2B_filefill7 = new double[f_lengthB[0]*(f_lengthB[2])];
    Real *x3B_filefill7 = new double[f_lengthB[0]*(f_lengthB[3])];
    if (should_fill_arr[7]) { 
        hdf5_open(fname_arr[7].c_str());
        hdf5_set_directory("/");
        hdf5_read_array(rho_filefill7, "prims.rho", 4, fill_fdims, fstart, fill_fdims, fill_fdims, fstart, H5T_IEEE_F64LE);
        hdf5_read_array(u_filefill7, "prims.u", 4, fill_fdims, fstart, fill_fdims, fill_fdims, fstart, H5T_IEEE_F64LE);
        hdf5_read_array(uvec_filefill7, "prims.uvec", 5, fill_fdims_vec, fstart_vec, fill_fdims_vec, fill_fdims_vec, fstart_vec, H5T_IEEE_F64LE);
        if (include_B) {
            if (b_ct) hdf5_read_array(B_filefill7, "cons.fB", 5, fill_fdimsB_vec, fstart_vec, fill_fdimsB_vec, fill_fdimsB_vec, fstart_vec,H5T_IEEE_F64LE);
            else if (pkgs.count("B_FluxCT")) hdf5_read_array(B_filefill7, "cons.B", 5, fill_fdims_vec, fstart_vec, fill_fdims_vec, fill_fdims_vec, fstart_vec,H5T_IEEE_F64LE);
        }
        hdf5_read_array(x1_filefill7, "VolumeLocations/x", 2, fill_fdims_x1, fstart_x, fill_fdims_x1, fill_fdims_x1, fstart_x, H5T_IEEE_F64LE);
        hdf5_read_array(x2_filefill7, "VolumeLocations/y", 2, fill_fdims_x2, fstart_x, fill_fdims_x2, fill_fdims_x2, fstart_x, H5T_IEEE_F64LE);
        hdf5_read_array(x3_filefill7, "VolumeLocations/z", 2, fill_fdims_x3, fstart_x, fill_fdims_x3, fill_fdims_x3, fstart_x, H5T_IEEE_F64LE);
        if (b_ct) {
            hdf5_read_array(x1B_filefill7, "Locations/x", 2, fill_fdimsB_x1, fstart_x, fill_fdimsB_x1, fill_fdimsB_x1, fstart_x, H5T_IEEE_F64LE);
            hdf5_read_array(x2B_filefill7, "Locations/y", 2, fill_fdimsB_x2, fstart_x, fill_fdimsB_x2, fill_fdimsB_x2, fstart_x, H5T_IEEE_F64LE);
            hdf5_read_array(x3B_filefill7, "Locations/z", 2, fill_fdimsB_x3, fstart_x, fill_fdimsB_x3, fill_fdimsB_x3, fstart_x, H5T_IEEE_F64LE);
        }
        hdf5_close();
    }

    // save the grid coordinate values to host array
    for (int iblocktemp = 0; iblocktemp < length[0]; iblocktemp++) {
        for (int itemp = 0; itemp < length[1]; itemp++) {
            x1_f_host(iblocktemp,itemp) = x1_file[length[1]*iblocktemp+itemp];
        }
        for (int jtemp = 0; jtemp < length[2]; jtemp++) {
            x2_f_host(iblocktemp,jtemp) = x2_file[length[2]*iblocktemp+jtemp];
        }
        for (int ktemp = 0; ktemp < length[3]; ktemp++) {
            x3_f_host(iblocktemp,ktemp) = x3_file[length[3]*iblocktemp+ktemp];
        }
    }
    // re-arrange uvec such that it can be read in the VLOOP
    int vector_file_index, scalar_file_index;
    for (int iblocktemp = 0; iblocktemp < length[0]; iblocktemp++) {
        for (int itemp = 0; itemp < length[1]; itemp++) {
            for (int jtemp = 0; jtemp < length[2]; jtemp++) {
                for (int ktemp = 0; ktemp < length[3]; ktemp++) {
                    scalar_file_index = length[1]*(length[2]*(length[3]*iblocktemp+ktemp)+jtemp)+itemp;

                    rho_f_host(iblocktemp,ktemp,jtemp,itemp) = rho_file[scalar_file_index];
                    u_f_host(iblocktemp,ktemp,jtemp,itemp) = u_file[scalar_file_index];
                    for (int ltemp = 0; ltemp < 3; ltemp++) {
                        //vector_file_index = 3*(scalar_file_index)+ltemp; // outdated parthenon phdf5 saving order
                        vector_file_index = length[1]*(length[2]*(length[3]*(NVEC*iblocktemp+ltemp)+ktemp)+jtemp)+itemp;
                        
                        uvec_f_host(ltemp,iblocktemp,ktemp,jtemp,itemp) = uvec_file[vector_file_index];
                        //if (include_B) B_f_host(ltemp,iblocktemp,ktemp,jtemp,itemp) = B_file[vector_file_index];
                    }
                }
            }
        }
    }
    if (include_B) {
        for (int iblocktemp = 0; iblocktemp < lengthB[0]; iblocktemp++) {
            for (int itemp = 0; itemp < lengthB[1]; itemp++) {
                x1B_f_host(iblocktemp, itemp) = x1B_file[lengthB[1] * iblocktemp + itemp];
                for (int jtemp = 0; jtemp < lengthB[2]; jtemp++) {
                    for (int ktemp = 0; ktemp < lengthB[3]; ktemp++) {
                        for (int ltemp = 0; ltemp < 3; ltemp++) {
                            vector_file_index = lengthB[1] * (lengthB[2] * (lengthB[3] * (NVEC * iblocktemp + ltemp) + ktemp) + jtemp) + itemp;
                            B_f_host(ltemp, iblocktemp, ktemp, jtemp, itemp) = B_file[vector_file_index];
                        }
                    }
                }
            }
            for (int jtemp = 0; jtemp < lengthB[2]; jtemp++)
                x2B_f_host(iblocktemp, jtemp) = x2B_file[lengthB[2] * iblocktemp + jtemp];
            for (int ktemp = 0; ktemp < lengthB[3]; ktemp++) 
                x3B_f_host(iblocktemp, ktemp) = x3B_file[lengthB[3] * iblocktemp + ktemp];
        }
    }
    
    // now repeat for fillfiles
    // save the grid coordinate values to host array
    for (int iblocktemp = 0; iblocktemp < f_length[0]; iblocktemp++) {
        for (int itemp = 0; itemp < f_length[1]; itemp++) {
            if (should_fill_arr[1]) x1_fill1_host(iblocktemp,itemp) = x1_filefill1[f_length[1]*iblocktemp+itemp];
            if (should_fill_arr[2]) x1_fill2_host(iblocktemp,itemp) = x1_filefill2[f_length[1]*iblocktemp+itemp];
            if (should_fill_arr[3]) x1_fill3_host(iblocktemp,itemp) = x1_filefill3[f_length[1]*iblocktemp+itemp];
            if (should_fill_arr[4]) x1_fill4_host(iblocktemp,itemp) = x1_filefill4[f_length[1]*iblocktemp+itemp];
            if (should_fill_arr[5]) x1_fill5_host(iblocktemp,itemp) = x1_filefill5[f_length[1]*iblocktemp+itemp];
            if (should_fill_arr[6]) x1_fill6_host(iblocktemp,itemp) = x1_filefill6[f_length[1]*iblocktemp+itemp];
            if (should_fill_arr[7]) x1_fill7_host(iblocktemp,itemp) = x1_filefill7[f_length[1]*iblocktemp+itemp];
        } for (int jtemp = 0; jtemp < f_length[2]; jtemp++) {
            if (should_fill_arr[1]) x2_fill1_host(iblocktemp,jtemp) = x2_filefill1[f_length[2]*iblocktemp+jtemp];
            if (should_fill_arr[2]) x2_fill2_host(iblocktemp,jtemp) = x2_filefill2[f_length[2]*iblocktemp+jtemp];
            if (should_fill_arr[3]) x2_fill3_host(iblocktemp,jtemp) = x2_filefill3[f_length[2]*iblocktemp+jtemp];
            if (should_fill_arr[4]) x2_fill4_host(iblocktemp,jtemp) = x2_filefill4[f_length[2]*iblocktemp+jtemp];
            if (should_fill_arr[5]) x2_fill5_host(iblocktemp,jtemp) = x2_filefill5[f_length[2]*iblocktemp+jtemp];
            if (should_fill_arr[6]) x2_fill6_host(iblocktemp,jtemp) = x2_filefill6[f_length[2]*iblocktemp+jtemp];
            if (should_fill_arr[7]) x2_fill7_host(iblocktemp,jtemp) = x2_filefill7[f_length[2]*iblocktemp+jtemp];
        } for (int ktemp = 0; ktemp < f_length[3]; ktemp++) {
            if (should_fill_arr[1]) x3_fill1_host(iblocktemp,ktemp) = x3_filefill1[f_length[3]*iblocktemp+ktemp];
            if (should_fill_arr[2]) x3_fill2_host(iblocktemp,ktemp) = x3_filefill2[f_length[3]*iblocktemp+ktemp];
            if (should_fill_arr[3]) x3_fill3_host(iblocktemp,ktemp) = x3_filefill3[f_length[3]*iblocktemp+ktemp];
            if (should_fill_arr[4]) x3_fill4_host(iblocktemp,ktemp) = x3_filefill4[f_length[3]*iblocktemp+ktemp];
            if (should_fill_arr[5]) x3_fill5_host(iblocktemp,ktemp) = x3_filefill5[f_length[3]*iblocktemp+ktemp];
            if (should_fill_arr[6]) x3_fill6_host(iblocktemp,ktemp) = x3_filefill6[f_length[3]*iblocktemp+ktemp];
            if (should_fill_arr[7]) x3_fill7_host(iblocktemp,ktemp) = x3_filefill7[f_length[3]*iblocktemp+ktemp];
        }
    }
    
    // re-arrange uvec such that it can be read in the VLOOP
    for (int iblocktemp = 0; iblocktemp < f_length[0]; iblocktemp++) {
        for (int itemp = 0; itemp < f_length[1]; itemp++) {
            for (int jtemp = 0; jtemp < f_length[2]; jtemp++) {
                for (int ktemp = 0; ktemp < f_length[3]; ktemp++) {
                    scalar_file_index = f_length[1]*(f_length[2]*(f_length[3]*iblocktemp+ktemp)+jtemp)+itemp;
                    if (should_fill_arr[1]) {
                        rho_fill1_host(iblocktemp,ktemp,jtemp,itemp) = rho_filefill1[scalar_file_index];
                        u_fill1_host(iblocktemp,ktemp,jtemp,itemp) = u_filefill1[scalar_file_index];
                    }
                    if (should_fill_arr[2]) {
                        rho_fill2_host(iblocktemp,ktemp,jtemp,itemp) = rho_filefill2[scalar_file_index];
                        u_fill2_host(iblocktemp,ktemp,jtemp,itemp) = u_filefill2[scalar_file_index];
                    }
                    if (should_fill_arr[3]) {
                        rho_fill3_host(iblocktemp,ktemp,jtemp,itemp) = rho_filefill3[scalar_file_index];
                        u_fill3_host(iblocktemp,ktemp,jtemp,itemp) = u_filefill3[scalar_file_index];
                    }
                    if (should_fill_arr[4]) {
                        rho_fill4_host(iblocktemp,ktemp,jtemp,itemp) = rho_filefill4[scalar_file_index];
                        u_fill4_host(iblocktemp,ktemp,jtemp,itemp) = u_filefill4[scalar_file_index];
                    }
                    if (should_fill_arr[5]) {
                        rho_fill5_host(iblocktemp,ktemp,jtemp,itemp) = rho_filefill5[scalar_file_index];
                        u_fill5_host(iblocktemp,ktemp,jtemp,itemp) = u_filefill5[scalar_file_index];
                    }
                    if (should_fill_arr[6]) {
                        rho_fill6_host(iblocktemp,ktemp,jtemp,itemp) = rho_filefill6[scalar_file_index];
                        u_fill6_host(iblocktemp,ktemp,jtemp,itemp) = u_filefill6[scalar_file_index];
                    }
                    if (should_fill_arr[7]) {
                        rho_fill7_host(iblocktemp,ktemp,jtemp,itemp) = rho_filefill7[scalar_file_index];
                        u_fill7_host(iblocktemp,ktemp,jtemp,itemp) = u_filefill7[scalar_file_index];
                    }
                    for (int ltemp = 0; ltemp < 3; ltemp++) {
                        vector_file_index = f_length[1]*(f_length[2]*(f_length[3]*(3*iblocktemp+ltemp)+ktemp)+jtemp)+itemp;
                        if (should_fill_arr[1]) uvec_fill1_host(ltemp,iblocktemp,ktemp,jtemp,itemp) = uvec_filefill1[vector_file_index];
                        if (should_fill_arr[2]) uvec_fill2_host(ltemp,iblocktemp,ktemp,jtemp,itemp) = uvec_filefill2[vector_file_index];
                        if (should_fill_arr[3]) uvec_fill3_host(ltemp,iblocktemp,ktemp,jtemp,itemp) = uvec_filefill3[vector_file_index];
                        if (should_fill_arr[4]) uvec_fill4_host(ltemp,iblocktemp,ktemp,jtemp,itemp) = uvec_filefill4[vector_file_index];
                        if (should_fill_arr[5]) uvec_fill5_host(ltemp,iblocktemp,ktemp,jtemp,itemp) = uvec_filefill5[vector_file_index];
                        if (should_fill_arr[6]) uvec_fill6_host(ltemp,iblocktemp,ktemp,jtemp,itemp) = uvec_filefill6[vector_file_index];
                        if (should_fill_arr[7]) uvec_fill7_host(ltemp,iblocktemp,ktemp,jtemp,itemp) = uvec_filefill7[vector_file_index];
                    }
                }
            }
        }
    }
    if (include_B) {
        for (int iblocktemp = 0; iblocktemp < f_lengthB[0]; iblocktemp++) {
            for (int itemp = 0; itemp < f_lengthB[1]; itemp++) {
                scalar_file_index = f_lengthB[1] * iblocktemp + itemp;
                if (should_fill_arr[1]) x1B_fill1_host(iblocktemp, itemp) = x1B_filefill1[scalar_file_index];
                if (should_fill_arr[2]) x1B_fill2_host(iblocktemp, itemp) = x1B_filefill2[scalar_file_index];
                if (should_fill_arr[3]) x1B_fill3_host(iblocktemp, itemp) = x1B_filefill3[scalar_file_index];
                if (should_fill_arr[4]) x1B_fill4_host(iblocktemp, itemp) = x1B_filefill4[scalar_file_index];
                if (should_fill_arr[5]) x1B_fill5_host(iblocktemp, itemp) = x1B_filefill5[scalar_file_index];
                if (should_fill_arr[6]) x1B_fill6_host(iblocktemp, itemp) = x1B_filefill6[scalar_file_index];
                if (should_fill_arr[7]) x1B_fill7_host(iblocktemp, itemp) = x1B_filefill7[scalar_file_index];
                for (int jtemp = 0; jtemp < f_lengthB[2]; jtemp++) {
                    for (int ktemp = 0; ktemp < f_lengthB[3]; ktemp++) {
                        for (int ltemp = 0; ltemp < 3; ltemp++) {
                            vector_file_index = f_lengthB[1] * (f_lengthB[2] * (f_lengthB[3] * (NVEC * iblocktemp + ltemp) + ktemp) + jtemp) + itemp;
                            if (should_fill_arr[1]) B_fill1_host(ltemp, iblocktemp, ktemp, jtemp, itemp) = B_filefill1[vector_file_index];
                            if (should_fill_arr[2]) B_fill2_host(ltemp, iblocktemp, ktemp, jtemp, itemp) = B_filefill2[vector_file_index];
                            if (should_fill_arr[3]) B_fill3_host(ltemp, iblocktemp, ktemp, jtemp, itemp) = B_filefill3[vector_file_index];
                            if (should_fill_arr[4]) B_fill4_host(ltemp, iblocktemp, ktemp, jtemp, itemp) = B_filefill4[vector_file_index];
                            if (should_fill_arr[5]) B_fill5_host(ltemp, iblocktemp, ktemp, jtemp, itemp) = B_filefill5[vector_file_index];
                            if (should_fill_arr[6]) B_fill6_host(ltemp, iblocktemp, ktemp, jtemp, itemp) = B_filefill6[vector_file_index];
                            if (should_fill_arr[7]) B_fill7_host(ltemp, iblocktemp, ktemp, jtemp, itemp) = B_filefill7[vector_file_index];
                        }
                    }
                }
            }
            for (int jtemp = 0; jtemp < f_lengthB[2]; jtemp++) {
                scalar_file_index = f_lengthB[2] * iblocktemp + jtemp;
                if (should_fill_arr[1]) x2B_fill1_host(iblocktemp, jtemp) = x2B_filefill1[scalar_file_index];
                if (should_fill_arr[2]) x2B_fill2_host(iblocktemp, jtemp) = x2B_filefill2[scalar_file_index];
                if (should_fill_arr[3]) x2B_fill3_host(iblocktemp, jtemp) = x2B_filefill3[scalar_file_index];
                if (should_fill_arr[4]) x2B_fill4_host(iblocktemp, jtemp) = x2B_filefill4[scalar_file_index];
                if (should_fill_arr[5]) x2B_fill5_host(iblocktemp, jtemp) = x2B_filefill5[scalar_file_index];
                if (should_fill_arr[6]) x2B_fill6_host(iblocktemp, jtemp) = x2B_filefill6[scalar_file_index];
                if (should_fill_arr[7]) x2B_fill7_host(iblocktemp, jtemp) = x2B_filefill7[scalar_file_index];
            }
            for (int ktemp = 0; ktemp < f_lengthB[3]; ktemp++) {
                scalar_file_index = f_lengthB[3] * iblocktemp + ktemp;
                if (should_fill_arr[1]) x3B_fill1_host(iblocktemp, ktemp) = x3B_filefill1[scalar_file_index];
                if (should_fill_arr[2]) x3B_fill2_host(iblocktemp, ktemp) = x3B_filefill2[scalar_file_index];
                if (should_fill_arr[3]) x3B_fill3_host(iblocktemp, ktemp) = x3B_filefill3[scalar_file_index];
                if (should_fill_arr[4]) x3B_fill4_host(iblocktemp, ktemp) = x3B_filefill4[scalar_file_index];
                if (should_fill_arr[5]) x3B_fill5_host(iblocktemp, ktemp) = x3B_filefill5[scalar_file_index];
                if (should_fill_arr[6]) x3B_fill6_host(iblocktemp, ktemp) = x3B_filefill6[scalar_file_index];
                if (should_fill_arr[7]) x3B_fill7_host(iblocktemp, ktemp) = x3B_filefill7[scalar_file_index];
            }
        }
    }

    // Deep copy to device
    x1_f_device.DeepCopy(x1_f_host);
    x2_f_device.DeepCopy(x2_f_host);
    x3_f_device.DeepCopy(x3_f_host);
    rho_f_device.DeepCopy(rho_f_host);
    u_f_device.DeepCopy(u_f_host);
    uvec_f_device.DeepCopy(uvec_f_host);
    if (include_B) {
        B_f_device.DeepCopy(B_f_host);
        if (b_ct) {
            x1B_f_device.DeepCopy(x1B_f_host);
            x2B_f_device.DeepCopy(x2B_f_host);
            x3B_f_device.DeepCopy(x3B_f_host);
        } else { // otherwise, same as cell-centered values
            x1B_f_device.DeepCopy(x1_f_host);
            x2B_f_device.DeepCopy(x2_f_host);
            x3B_f_device.DeepCopy(x3_f_host);
        }
    }
    if (should_fill_arr[1]) {
        x1_fill1_device.DeepCopy(x1_fill1_host);
        x2_fill1_device.DeepCopy(x2_fill1_host);
        x3_fill1_device.DeepCopy(x3_fill1_host);
        rho_fill1_device.DeepCopy(rho_fill1_host);
        u_fill1_device.DeepCopy(u_fill1_host);
        uvec_fill1_device.DeepCopy(uvec_fill1_host);
        if (include_B) {
            B_fill1_device.DeepCopy(B_fill1_host);
            if (b_ct) {
                x1B_fill1_device.DeepCopy(x1B_fill1_host);
                x2B_fill1_device.DeepCopy(x2B_fill1_host);
                x3B_fill1_device.DeepCopy(x3B_fill1_host);
            } else {
                x1B_fill1_device.DeepCopy(x1_fill1_host);
                x2B_fill1_device.DeepCopy(x2_fill1_host);
                x3B_fill1_device.DeepCopy(x3_fill1_host);
            }
        }
    }
    if (should_fill_arr[2]) {
        x1_fill2_device.DeepCopy(x1_fill2_host);
        x2_fill2_device.DeepCopy(x2_fill2_host);
        x3_fill2_device.DeepCopy(x3_fill2_host);
        rho_fill2_device.DeepCopy(rho_fill2_host);
        u_fill2_device.DeepCopy(u_fill2_host);
        uvec_fill2_device.DeepCopy(uvec_fill2_host);
        if (include_B) {
            B_fill2_device.DeepCopy(B_fill2_host);
            if (b_ct) {
                x1B_fill2_device.DeepCopy(x1B_fill2_host);
                x2B_fill2_device.DeepCopy(x2B_fill2_host);
                x3B_fill2_device.DeepCopy(x3B_fill2_host);
            } else {
                x1B_fill2_device.DeepCopy(x1_fill2_host);
                x2B_fill2_device.DeepCopy(x2_fill2_host);
                x3B_fill2_device.DeepCopy(x3_fill2_host);
            }
        }
    }
    if (should_fill_arr[3]) {
        x1_fill3_device.DeepCopy(x1_fill3_host);
        x2_fill3_device.DeepCopy(x2_fill3_host);
        x3_fill3_device.DeepCopy(x3_fill3_host);
        rho_fill3_device.DeepCopy(rho_fill3_host);
        u_fill3_device.DeepCopy(u_fill3_host);
        uvec_fill3_device.DeepCopy(uvec_fill3_host);
        if (include_B) {
            B_fill3_device.DeepCopy(B_fill3_host);
            if (b_ct) {
                x1B_fill3_device.DeepCopy(x1B_fill3_host);
                x2B_fill3_device.DeepCopy(x2B_fill3_host);
                x3B_fill3_device.DeepCopy(x3B_fill3_host);
            } else {
                x1B_fill3_device.DeepCopy(x1_fill3_host);
                x2B_fill3_device.DeepCopy(x2_fill3_host);
                x3B_fill3_device.DeepCopy(x3_fill3_host);
            }
        }
    }
    if (should_fill_arr[4]) {
        x1_fill4_device.DeepCopy(x1_fill4_host);
        x2_fill4_device.DeepCopy(x2_fill4_host);
        x3_fill4_device.DeepCopy(x3_fill4_host);
        rho_fill4_device.DeepCopy(rho_fill4_host);
        u_fill4_device.DeepCopy(u_fill4_host);
        uvec_fill4_device.DeepCopy(uvec_fill4_host);
        if (include_B) {
            B_fill4_device.DeepCopy(B_fill4_host);
            if (b_ct) {
                x1B_fill4_device.DeepCopy(x1B_fill4_host);
                x2B_fill4_device.DeepCopy(x2B_fill4_host);
                x3B_fill4_device.DeepCopy(x3B_fill4_host);
            } else {
                x1B_fill4_device.DeepCopy(x1_fill4_host);
                x2B_fill4_device.DeepCopy(x2_fill4_host);
                x3B_fill4_device.DeepCopy(x3_fill4_host);
            }
        }
    }
    if (should_fill_arr[5]) {
        x1_fill5_device.DeepCopy(x1_fill5_host);
        x2_fill5_device.DeepCopy(x2_fill5_host);
        x3_fill5_device.DeepCopy(x3_fill5_host);
        rho_fill5_device.DeepCopy(rho_fill5_host);
        u_fill5_device.DeepCopy(u_fill5_host);
        uvec_fill5_device.DeepCopy(uvec_fill5_host);
        if (include_B) {
            B_fill5_device.DeepCopy(B_fill5_host);
            if (b_ct) {
                x1B_fill5_device.DeepCopy(x1B_fill5_host);
                x2B_fill5_device.DeepCopy(x2B_fill5_host);
                x3B_fill5_device.DeepCopy(x3B_fill5_host);
            } else {
                x1B_fill5_device.DeepCopy(x1_fill5_host);
                x2B_fill5_device.DeepCopy(x2_fill5_host);
                x3B_fill5_device.DeepCopy(x3_fill5_host);
            }
        }
    }
    if (should_fill_arr[6]) {
        x1_fill6_device.DeepCopy(x1_fill6_host);
        x2_fill6_device.DeepCopy(x2_fill6_host);
        x3_fill6_device.DeepCopy(x3_fill6_host);
        rho_fill6_device.DeepCopy(rho_fill6_host);
        u_fill6_device.DeepCopy(u_fill6_host);
        uvec_fill6_device.DeepCopy(uvec_fill6_host);
        if (include_B) {
            B_fill6_device.DeepCopy(B_fill6_host);
            if (b_ct) {
                x1B_fill6_device.DeepCopy(x1B_fill6_host);
                x2B_fill6_device.DeepCopy(x2B_fill6_host);
                x3B_fill6_device.DeepCopy(x3B_fill6_host);
            } else {
                x1B_fill6_device.DeepCopy(x1_fill6_host);
                x2B_fill6_device.DeepCopy(x2_fill6_host);
                x3B_fill6_device.DeepCopy(x3_fill6_host);
            }
        }
    }
    if (should_fill_arr[7]) {
        x1_fill7_device.DeepCopy(x1_fill7_host);
        x2_fill7_device.DeepCopy(x2_fill7_host);
        x3_fill7_device.DeepCopy(x3_fill7_host);
        rho_fill7_device.DeepCopy(rho_fill7_host);
        u_fill7_device.DeepCopy(u_fill7_host);
        uvec_fill7_device.DeepCopy(uvec_fill7_host);
        if (include_B) {
            B_fill7_device.DeepCopy(B_fill7_host);
            if (b_ct) {
                x1B_fill7_device.DeepCopy(x1B_fill7_host);
                x2B_fill7_device.DeepCopy(x2B_fill7_host);
                x3B_fill7_device.DeepCopy(x3B_fill7_host);
            } else {
                x1B_fill7_device.DeepCopy(x1_fill7_host);
                x2B_fill7_device.DeepCopy(x2_fill7_host);
                x3B_fill7_device.DeepCopy(x3_fill7_host);
            }
        }
    }
    
    Kokkos::fence();

    PackIndexMap prims_map;
    auto P = GRMHD::PackMHDPrims(rc.get(), prims_map);
    const VarMap m_p(prims_map, false);

    // Device-side interpolate & copy into the mirror array
    if (MPIRank0() && verbose > 0) {
        std::cout << "Initializing KHARMA restart.  Filling from " << fname_arr[0]
                    << ", " << fname_arr[1] << ", ..." << std::endl;
    }

    // Read to the entire meshblock -- we'll set the Dirichlet boundaries based on the
    // ghost zone data we read here.
    auto domain = IndexDomain::entire;
    int is = pmb->cellbounds.is(domain), ie = pmb->cellbounds.ie(domain);
    int js = pmb->cellbounds.js(domain), je = pmb->cellbounds.je(domain);
    int ks = pmb->cellbounds.ks(domain), ke = pmb->cellbounds.ke(domain);

    pmb->par_for("copy_restart_state_kharma", ks, ke, js, je, is, ie,
        KOKKOS_LAMBDA (const int &k, const int &j, const int &i) {
            GReal X[GR_DIM];
            G.coord(k, j, i, Loci::center, X);
            if ((X[1] >= fx1min_arr[0] - fnghost * dx[1]) && (X[1] <=  fx1max_arr[0] + fnghost * dx[1])) { // fill with the fname
                get_prim_restart_kharma(X, P, m_p, dx, length,
                    x1_f_device, x2_f_device, x3_f_device, rho_f_device, u_f_device, uvec_f_device,
                    k, j, i);
            } else if ((should_fill_arr[1]) && (X[1] >= fx1min_arr[1] - fnghost * dx[1]) && (X[1] <= fx1max_arr[1] + fnghost * dx[1])) { // fill with the fname_fill1
                get_prim_restart_kharma(X, P, m_p, dx, f_length,
                    x1_fill1_device, x2_fill1_device, x3_fill1_device, rho_fill1_device, u_fill1_device, uvec_fill1_device,
                    k, j, i);
            } else if ((should_fill_arr[2]) && (X[1] >= fx1min_arr[2] - fnghost * dx[1]) && (X[1] <= fx1max_arr[2] + fnghost * dx[1])) { // fill with the fname_fill2
                get_prim_restart_kharma(X, P, m_p, dx, f_length,
                    x1_fill2_device, x2_fill2_device, x3_fill2_device, rho_fill2_device, u_fill2_device, uvec_fill2_device,
                    k, j, i);
            } else if ((should_fill_arr[3]) && (X[1] >= fx1min_arr[3] - fnghost * dx[1]) && (X[1] <= fx1max_arr[3] + fnghost * dx[1])) { // fill with the fname_fill3
                get_prim_restart_kharma(X, P, m_p, dx, f_length,
                    x1_fill3_device, x2_fill3_device, x3_fill3_device, rho_fill3_device, u_fill3_device, uvec_fill3_device,
                    k, j, i);
            } else if ((should_fill_arr[4]) && (X[1] >= fx1min_arr[4] - fnghost * dx[1]) && (X[1] <= fx1max_arr[4] + fnghost * dx[1])) { // fill with the fname_fill4
                get_prim_restart_kharma(X, P, m_p, dx, f_length,
                    x1_fill4_device, x2_fill4_device, x3_fill4_device, rho_fill4_device, u_fill4_device, uvec_fill4_device,
                    k, j, i);
            } else if ((should_fill_arr[5]) && (X[1] >= fx1min_arr[5] - fnghost * dx[1]) && (X[1] <= fx1max_arr[5] + fnghost * dx[1])) { // fill with the fname_fill5
                get_prim_restart_kharma(X, P, m_p, dx, f_length,
                    x1_fill5_device, x2_fill5_device, x3_fill5_device, rho_fill5_device, u_fill5_device, uvec_fill5_device,
                    k, j, i);
            } else if ((should_fill_arr[6]) && (X[1] >= fx1min_arr[6] - fnghost * dx[1]) && (X[1] <= fx1max_arr[6] + fnghost * dx[1])) { // fill with the fname_fill6
                get_prim_restart_kharma(X, P, m_p, dx, f_length,
                    x1_fill6_device, x2_fill6_device, x3_fill6_device, rho_fill6_device, u_fill6_device, uvec_fill6_device,
                    k, j, i);
            } else if (should_fill_arr[7]) { // fill with the fname_fill7
                get_prim_restart_kharma(X, P, m_p, dx, f_length,
                    x1_fill7_device, x2_fill7_device, x3_fill7_device, rho_fill7_device, u_fill7_device, uvec_fill7_device,
                    k, j, i);
            } else {
                printf("HYERIN: no corresponding file found to fill!\n");
            }
        }
    );
    if (include_B) {
        Loci loc;
        VLOOP {
            if (b_ct) loc = (Loci) v;
            else loc = Loci::center;
            pmb->par_for("copy_restart_state_kharma", ks, ke + (dir_of(loc) == 3), js, je + (dir_of(loc) == 2), is, ie + (dir_of(loc) == 1),
                KOKKOS_LAMBDA (const int &k, const int &j, const int &i) {
                    GReal X[GR_DIM];
                    G.coord(k, j, i, loc, X);
                    if ((X[1] >= fx1min_arr[0] - fnghost * dx[1]) && (X[1] <=  fx1max_arr[0] + fnghost * dx[1])) { // fill with the fname
                        get_B_restart_kharma(X, dx, loc, length, lengthB, 
                                x1_f_device, x2_f_device, x3_f_device, 
                                x1B_f_device, x2B_f_device, x3B_f_device, B_f_device, B_Save, v, k, j, i);
                    } else if ((should_fill_arr[1]) && (X[1] >= fx1min_arr[1] - fnghost * dx[1]) && (X[1] <= fx1max_arr[1] + fnghost * dx[1])) { // fill with the fname_fill1
                        get_B_restart_kharma(X, dx, loc, f_length, f_lengthB, 
                                x1_fill1_device, x2_fill1_device, x3_fill1_device, 
                                x1B_fill1_device, x2B_fill1_device, x3B_fill1_device, B_fill1_device, B_Save, v, k, j, i);
                    } else if ((should_fill_arr[2]) && (X[1] >= fx1min_arr[2] - fnghost * dx[1]) && (X[1] <= fx1max_arr[2] + fnghost * dx[1])) { // fill with the fname_fill2
                        get_B_restart_kharma(X, dx, loc, f_length, f_lengthB, 
                                x1_fill2_device, x2_fill2_device, x3_fill2_device, 
                                x1B_fill2_device, x2B_fill2_device, x3B_fill2_device, B_fill2_device, B_Save, v, k, j, i);
                    } else if ((should_fill_arr[3]) && (X[1] >= fx1min_arr[3] - fnghost * dx[1]) && (X[1] <= fx1max_arr[3] + fnghost * dx[1])) { // fill with the fname_fill3
                        get_B_restart_kharma(X, dx, loc, f_length, f_lengthB, 
                                x1_fill3_device, x2_fill3_device, x3_fill3_device, 
                                x1B_fill3_device, x2B_fill3_device, x3B_fill3_device, B_fill3_device, B_Save, v, k, j, i);
                    } else if ((should_fill_arr[4]) && (X[1] >= fx1min_arr[4] - fnghost * dx[1]) && (X[1] <= fx1max_arr[4] + fnghost * dx[1])) { // fill with the fname_fill4
                        get_B_restart_kharma(X, dx, loc, f_length, f_lengthB, 
                                x1_fill4_device, x2_fill4_device, x3_fill4_device, 
                                x1B_fill4_device, x2B_fill4_device, x3B_fill4_device, B_fill4_device, B_Save, v, k, j, i);
                    } else if ((should_fill_arr[5]) && (X[1] >= fx1min_arr[5] - fnghost * dx[1]) && (X[1] <= fx1max_arr[5] + fnghost * dx[1])) { // fill with the fname_fill5
                        get_B_restart_kharma(X, dx, loc, f_length, f_lengthB, 
                                x1_fill5_device, x2_fill5_device, x3_fill5_device, 
                                x1B_fill5_device, x2B_fill5_device, x3B_fill5_device, B_fill5_device, B_Save, v, k, j, i);
                    } else if ((should_fill_arr[6]) && (X[1] >= fx1min_arr[6] - fnghost * dx[1]) && (X[1] <= fx1max_arr[6] + fnghost * dx[1])) { // fill with the fname_fill6
                        get_B_restart_kharma(X, dx, loc, f_length, f_lengthB, 
                                x1_fill6_device, x2_fill6_device, x3_fill6_device, 
                                x1B_fill6_device, x2B_fill6_device, x3B_fill6_device, B_fill6_device, B_Save, v, k, j, i);
                    } else if (should_fill_arr[7]) { // fill with the fname_fill7
                        get_B_restart_kharma(X, dx, loc, f_length, f_lengthB, 
                                x1_fill7_device, x2_fill7_device, x3_fill7_device, 
                                x1B_fill7_device, x2B_fill7_device, x3B_fill7_device, B_fill7_device, B_Save, v, k, j, i);
                    } else {
                        printf("HYERIN: no corresponding file found to fill for B!\n");
                    }
                }
            );
        }
    }

    return TaskStatus::complete;
}
