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

#include "hdf5_utils.h"
#include "types.hpp"

#include <sys/stat.h>
#include <ctype.h>

//using namespace Kokkos; // Hyerin: 10/07/22 comment this out, use par_for instead


// TODO
// Record & read:
// 1. startx/stopx/dx
// 2. coordinate name FMKS/MKS/etc
// 3. all coordinate params in play
// 4. Electron MODEL bitflag param
// 5. nprim for sanity check?
// 6. Indication of EMHD vs MHD

// TODO this code is very specific to spherical systems/boundares or entirely periodic boxes.
// No other boundaries/geometries are really supported.
//
// Reads in KHARMA restart file but at a different simulation size

void ReadFillFile(int i, std::unique_ptr<ParameterInput>& pin) {
    char str[20];
    sprintf(str, "fname_fill%d", i);
    auto fname_fill = pin->GetOrAddString("resize_restart", str, "none");

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

        sprintf(str, "restart_x1min_f%d", i);
        pin->SetReal("parthenon/mesh", str, fx1min);

        sprintf(str, "restart_x1max_f%d", i);
        pin->SetReal("parthenon/mesh", str, fx1max);

        sprintf(str, "restart%d_nx1", i);
        pin->SetInteger("parthenon/mesh", str, fnx1);
        pin->SetInteger("parthenon/meshblock", str, fmbnx1);
        sprintf(str, "restart%d_nx2", i);
        pin->SetInteger("parthenon/mesh", str, fnx2);
        pin->SetInteger("parthenon/meshblock", str, fmbnx2);
        sprintf(str, "restart%d_nx3", i);
        pin->SetInteger("parthenon/mesh", str, fnx3);
        pin->SetInteger("parthenon/meshblock", str, fmbnx3);
    }
}

void ReadKharmaRestartHeader(std::string fname, std::unique_ptr<ParameterInput>& pin)
{
    bool use_dt = pin->GetOrAddBoolean("resize_restart", "use_dt", true);
    bool use_tf = pin->GetOrAddBoolean("resize_restart", "use_tf", false);

    // Read input from restart file 
    // (from external/parthenon/src/parthenon_manager.cpp)
    std::unique_ptr<RestartReader> restartReader;
    restartReader = std::make_unique<RestartReader>(fname.c_str());

    // Load input stream
    std::unique_ptr<ParameterInput> fpinput;
    fpinput = std::make_unique<ParameterInput>();
    auto inputString = restartReader->GetAttr<std::string>("Input", "File");
    std::istringstream is(inputString);
    fpinput->LoadFromStream(is);

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
    //pin->SetString("b_field", "type", fBfield); // (12/07/22) Hyerin need to test

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

    Real  a, hslope;//, Rout;
    a = fpinput->GetReal("coordinates", "a");
    pin->SetReal("coordinates", "a", a);
    hslope = fpinput->GetReal("coordinates", "hslope");
    pin->SetReal("coordinates", "hslope", hslope);

    // close hdf5 file to prevent HDF5 hangs and corrupted files
    restartReader = nullptr;


    ReadFillFile(1, pin);
    ReadFillFile(2, pin);
    ReadFillFile(3, pin);
    ReadFillFile(4, pin);
    ReadFillFile(5, pin);
    ReadFillFile(6, pin);
    ReadFillFile(7, pin);
    /*
    Real fx1min_f1 = 0.;
    Real fx1max_f1 = 0.;
    auto fname_fill1 = pin->GetOrAddString("resize_restart", "fname_fill1", "none");
    if (!(fname_fill1 == "none")) {
        std::unique_ptr<RestartReader> restartReader_f1;
        restartReader_f1 = std::make_unique<RestartReader>(fname_fill1.c_str());

        // Load input stream
        std::unique_ptr<ParameterInput> fpinput_f1;
        fpinput_f1 = std::make_unique<ParameterInput>();
        auto inputString_f1 = restartReader_f1->GetAttr<std::string>("Input", "File");
        std::istringstream is_f1(inputString_f1);
        fpinput_f1->LoadFromStream(is_f1);

        fx1min_f1 = fpinput_f1->GetReal("parthenon/mesh", "x1min");
        fx1max_f1 = fpinput_f1->GetReal("parthenon/mesh", "x1max");
        
        restartReader_f1 = nullptr;
    }
    pin->SetReal("parthenon/mesh", "restart_x1min_f1", fx1min_f1);
    pin->SetReal("parthenon/mesh", "restart_x1max_f1", fx1max_f1);
    */
}

TaskStatus ReadKharmaRestart(MeshBlockData<Real> *rc, ParameterInput *pin)
{
    Flag(rc, "Restarting from KHARMA checkpoint file");

    auto pmb = rc->GetBlockPointer();

    const int n1tot = pin->GetInteger("parthenon/mesh", "restart_nx1");
    const int n2tot = pin->GetInteger("parthenon/mesh", "restart_nx2");
    const int n3tot = pin->GetInteger("parthenon/mesh", "restart_nx3");
    const int n1mb = pin->GetInteger("parthenon/meshblock", "restart_nx1");
    const int n2mb = pin->GetInteger("parthenon/meshblock", "restart_nx2");
    const int n3mb = pin->GetInteger("parthenon/meshblock", "restart_nx3");
    auto fname = pin->GetString("resize_restart", "fname"); // Require this, don't guess. 1st priority
    auto fname_fill1 = pin->GetOrAddString("resize_restart", "fname_fill1", "none"); // 2nd
    auto fname_fill2 = pin->GetOrAddString("resize_restart", "fname_fill2", "none"); // 3rd
    auto fname_fill3 = pin->GetOrAddString("resize_restart", "fname_fill3", "none"); // 4th
    auto fname_fill4 = pin->GetOrAddString("resize_restart", "fname_fill4", "none"); // 5th
    auto fname_fill5 = pin->GetOrAddString("resize_restart", "fname_fill5", "none"); // 6th
    auto fname_fill6 = pin->GetOrAddString("resize_restart", "fname_fill6", "none"); // 7th
    auto fname_fill7 = pin->GetOrAddString("resize_restart", "fname_fill7", "none"); // 8th
    const Real fx1min = pin->GetReal("parthenon/mesh", "restart_x1min");
    const Real fx1max = pin->GetReal("parthenon/mesh", "restart_x1max");
    const Real fx1min_f1 = pin->GetOrAddReal("parthenon/mesh", "restart_x1min_f1", -1);
    const Real fx1max_f1 = pin->GetOrAddReal("parthenon/mesh", "restart_x1max_f1", -1);
    const Real fx1min_f2 = pin->GetOrAddReal("parthenon/mesh", "restart_x1min_f2", -1);
    const Real fx1max_f2 = pin->GetOrAddReal("parthenon/mesh", "restart_x1max_f2", -1);
    const Real fx1min_f3 = pin->GetOrAddReal("parthenon/mesh", "restart_x1min_f3", -1);
    const Real fx1max_f3 = pin->GetOrAddReal("parthenon/mesh", "restart_x1max_f3", -1);
    const Real fx1min_f4 = pin->GetOrAddReal("parthenon/mesh", "restart_x1min_f4", -1);
    const Real fx1max_f4 = pin->GetOrAddReal("parthenon/mesh", "restart_x1max_f4", -1);
    const Real fx1min_f5 = pin->GetOrAddReal("parthenon/mesh", "restart_x1min_f5", -1);
    const Real fx1max_f5 = pin->GetOrAddReal("parthenon/mesh", "restart_x1max_f5", -1);
    const Real fx1min_f6 = pin->GetOrAddReal("parthenon/mesh", "restart_x1min_f6", -1);
    const Real fx1max_f6 = pin->GetOrAddReal("parthenon/mesh", "restart_x1max_f6", -1);
    const Real fx1min_f7 = pin->GetOrAddReal("parthenon/mesh", "restart_x1min_f7", -1);
    const Real fx1max_f7 = pin->GetOrAddReal("parthenon/mesh", "restart_x1max_f7", -1);
    // TODO Hyerin (07/23/23) assume all the "fill" files have the same dimension.
    const int f_n1tot = pin->GetOrAddInteger("parthenon/mesh", "restart1_nx1", -1);
    const int f_n2tot = pin->GetOrAddInteger("parthenon/mesh", "restart1_nx2", -1);
    const int f_n3tot = pin->GetOrAddInteger("parthenon/mesh", "restart1_nx3", -1);
    const int f_n1mb = pin->GetOrAddInteger("parthenon/meshblock", "restart1_nx1", -1);
    const int f_n2mb = pin->GetOrAddInteger("parthenon/meshblock", "restart1_nx2", -1);
    const int f_n3mb = pin->GetOrAddInteger("parthenon/meshblock", "restart1_nx3", -1);
    const Real x1min = pin->GetReal("parthenon/mesh", "x1min");
    const Real x2min = pin->GetReal("parthenon/mesh", "x2min");
    const Real x2max = pin->GetReal("parthenon/mesh", "x2max");
    const Real x3min = pin->GetReal("parthenon/mesh", "x3min");
    const Real x3max = pin->GetReal("parthenon/mesh", "x3max");
    const int nghost = pin->GetReal("parthenon/mesh", "restart_nghost");
    const bool ghost_zones = pin->GetBoolean("parthenon/mesh", "restart_ghostzones");
    auto fBfield = pin->GetOrAddString("b_field", "type", "none");
    const Real uphi = pin->GetOrAddReal("bondi", "uphi", 0.); 
    const Real vacuum_logrho = pin->GetOrAddReal("bondi", "vacuum_logrho", 0.);
    const Real vacuum_log_u_over_rho = pin->GetOrAddReal("bondi", "vacuum_log_u_over_rho", 0.);

    // Add these to package properties, since they continue to be needed on boundaries
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("rnx1")))
        pmb->packages.Get("GRMHD")->AddParam<int>("rnx1", n1tot);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("rnx2")))
        pmb->packages.Get("GRMHD")->AddParam<int>("rnx2", n2tot);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("rnx3")))
        pmb->packages.Get("GRMHD")->AddParam<int>("rnx3", n3tot);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("rmbnx1")))
        pmb->packages.Get("GRMHD")->AddParam<int>("rmbnx1", n1mb);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("rmbnx2")))
        pmb->packages.Get("GRMHD")->AddParam<int>("rmbnx2", n2mb);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("rmbnx3")))
        pmb->packages.Get("GRMHD")->AddParam<int>("rmbnx3", n3mb);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("rfnx1")))
        pmb->packages.Get("GRMHD")->AddParam<int>("rfnx1", f_n1tot);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("rfnx2")))
        pmb->packages.Get("GRMHD")->AddParam<int>("rfnx2", f_n2tot);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("rfnx3")))
        pmb->packages.Get("GRMHD")->AddParam<int>("rfnx3", f_n3tot);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("rfmbnx1")))
        pmb->packages.Get("GRMHD")->AddParam<int>("rfmbnx1", f_n1mb);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("rfmbnx2")))
        pmb->packages.Get("GRMHD")->AddParam<int>("rfmbnx2", f_n2mb);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("rfmbnx3")))
        pmb->packages.Get("GRMHD")->AddParam<int>("rfmbnx3", f_n3mb);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("fname")))
        pmb->packages.Get("GRMHD")->AddParam<std::string>("fname", fname);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("fname_fill1")))
        pmb->packages.Get("GRMHD")->AddParam<std::string>("fname_fill1", fname_fill1);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("fname_fill2")))
        pmb->packages.Get("GRMHD")->AddParam<std::string>("fname_fill2", fname_fill2);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("fname_fill3")))
        pmb->packages.Get("GRMHD")->AddParam<std::string>("fname_fill3", fname_fill3);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("fname_fill4")))
        pmb->packages.Get("GRMHD")->AddParam<std::string>("fname_fill4", fname_fill4);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("fname_fill5")))
        pmb->packages.Get("GRMHD")->AddParam<std::string>("fname_fill5", fname_fill5);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("fname_fill6")))
        pmb->packages.Get("GRMHD")->AddParam<std::string>("fname_fill6", fname_fill6);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("fname_fill7")))
        pmb->packages.Get("GRMHD")->AddParam<std::string>("fname_fill7", fname_fill7);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("rx1min")))
        pmb->packages.Get("GRMHD")->AddParam<Real>("rx1min", fx1min);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("rx1max")))
        pmb->packages.Get("GRMHD")->AddParam<Real>("rx1max", fx1max);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("rx1min_f1")))
        pmb->packages.Get("GRMHD")->AddParam<Real>("rx1min_f1", fx1min_f1);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("rx1max_f1")))
        pmb->packages.Get("GRMHD")->AddParam<Real>("rx1max_f1", fx1max_f1);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("rx1min_f2")))
        pmb->packages.Get("GRMHD")->AddParam<Real>("rx1min_f2", fx1min_f2);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("rx1max_f2")))
        pmb->packages.Get("GRMHD")->AddParam<Real>("rx1max_f2", fx1max_f2);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("rx1min_f3")))
        pmb->packages.Get("GRMHD")->AddParam<Real>("rx1min_f3", fx1min_f3);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("rx1max_f3")))
        pmb->packages.Get("GRMHD")->AddParam<Real>("rx1max_f3", fx1max_f3);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("rx1min_f4")))
        pmb->packages.Get("GRMHD")->AddParam<Real>("rx1min_f4", fx1min_f4);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("rx1max_f4")))
        pmb->packages.Get("GRMHD")->AddParam<Real>("rx1max_f4", fx1max_f4);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("rx1min_f5")))
        pmb->packages.Get("GRMHD")->AddParam<Real>("rx1min_f5", fx1min_f5);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("rx1max_f5")))
        pmb->packages.Get("GRMHD")->AddParam<Real>("rx1max_f5", fx1max_f5);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("rx1min_f6")))
        pmb->packages.Get("GRMHD")->AddParam<Real>("rx1min_f6", fx1min_f6);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("rx1max_f6")))
        pmb->packages.Get("GRMHD")->AddParam<Real>("rx1max_f6", fx1max_f6);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("rx1min_f7")))
        pmb->packages.Get("GRMHD")->AddParam<Real>("rx1min_f7", fx1min_f7);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("rx1max_f7")))
        pmb->packages.Get("GRMHD")->AddParam<Real>("rx1max_f7", fx1max_f7);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("x1min")))
        pmb->packages.Get("GRMHD")->AddParam<Real>("x1min", x1min);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("x2min")))
        pmb->packages.Get("GRMHD")->AddParam<Real>("x2min", x2min);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("x2max")))
        pmb->packages.Get("GRMHD")->AddParam<Real>("x2max", x2max);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("x3min")))
        pmb->packages.Get("GRMHD")->AddParam<Real>("x3min", x3min);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("x3max")))
        pmb->packages.Get("GRMHD")->AddParam<Real>("x3max", x3max);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("rnghost")))
        pmb->packages.Get("GRMHD")->AddParam<int>("rnghost", nghost);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("rghostzones")))
        pmb->packages.Get("GRMHD")->AddParam<bool>("rghostzones", ghost_zones);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("b_field_type")))
        pmb->packages.Get("GRMHD")->AddParam<std::string>("b_field_type", fBfield);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("uphi")))
        pmb->packages.Get("GRMHD")->AddParam<Real>("uphi", uphi);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("vacuum_logrho")))
        pmb->packages.Get("GRMHD")->AddParam<Real>("vacuum_logrho", vacuum_logrho);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("vacuum_log_u_over_rho")))
        pmb->packages.Get("GRMHD")->AddParam<Real>("vacuum_log_u_over_rho", vacuum_log_u_over_rho);

    // Set the whole domain
    SetKharmaRestart(rc);

   return TaskStatus::complete;
}

TaskStatus SetKharmaRestart(MeshBlockData<Real> *rc, IndexDomain domain, bool coarse)
{
    Flag(rc, "Setting KHARMA restart zones");
    auto pmb = rc->GetBlockPointer();
    auto b_field_type = pmb->packages.Get("GRMHD")->Param<std::string>("b_field_type");
    const bool include_B = (b_field_type != "none");
    // A placeholder to save the B fields for SeedBField
    GridVector B_Save;
    if (include_B) B_Save = rc->Get("B_Save").data;
    const Real uphi = pmb->packages.Get("GRMHD")->Param<Real>("uphi");
    const Real vacuum_logrho = pmb->packages.Get("GRMHD")->Param<Real>("vacuum_logrho");
    const Real vacuum_log_u_over_rho = pmb->packages.Get("GRMHD")->Param<Real>("vacuum_log_u_over_rho");

    auto& G = pmb->coords;
    
    // Size/domain of the MeshBlock we're reading to
    int is, ie;
    if (domain == IndexDomain::outer_x1) {// copying from bondi
        is = pmb->cellbounds.GetBoundsI(IndexDomain::interior).e+1;
        ie = pmb->cellbounds.GetBoundsI(IndexDomain::entire).e;
    } else if (domain == IndexDomain::inner_x1) {
        is = pmb->cellbounds.GetBoundsI(IndexDomain::entire).s;
        ie = pmb->cellbounds.GetBoundsI(IndexDomain::interior).s-1;
    } else {
        is = pmb->cellbounds.is(domain);
        ie = pmb->cellbounds.ie(domain);
    }
    int js = pmb->cellbounds.js(domain), je = pmb->cellbounds.je(domain);
    int ks = pmb->cellbounds.ks(domain), ke = pmb->cellbounds.ke(domain);
    //IndexRange block = IndexRange{0, nb - 1};
    
    const int n1tot = pmb->packages.Get("GRMHD")->Param<int>("rnx1");
    const int n2tot = pmb->packages.Get("GRMHD")->Param<int>("rnx2");
    const int n3tot = pmb->packages.Get("GRMHD")->Param<int>("rnx3");
    hsize_t n1mb = pmb->packages.Get("GRMHD")->Param<int>("rmbnx1");
    hsize_t n2mb = pmb->packages.Get("GRMHD")->Param<int>("rmbnx2");
    hsize_t n3mb = pmb->packages.Get("GRMHD")->Param<int>("rmbnx3");
    hsize_t nBlocks = (int) (n1tot*n2tot*n3tot)/(n1mb*n2mb*n3mb);
    const int f_n1tot = pmb->packages.Get("GRMHD")->Param<int>("rfnx1");
    const int f_n2tot = pmb->packages.Get("GRMHD")->Param<int>("rfnx2");
    const int f_n3tot = pmb->packages.Get("GRMHD")->Param<int>("rfnx3");
    hsize_t f_n1mb = pmb->packages.Get("GRMHD")->Param<int>("rfmbnx1");
    hsize_t f_n2mb = pmb->packages.Get("GRMHD")->Param<int>("rfmbnx2");
    hsize_t f_n3mb = pmb->packages.Get("GRMHD")->Param<int>("rfmbnx3");
    hsize_t f_nBlocks = (int) (f_n1tot*f_n2tot*f_n3tot)/(f_n1mb*f_n2mb*f_n3mb);
    auto fname = pmb->packages.Get("GRMHD")->Param<std::string>("fname");
    auto fname_fill1 = pmb->packages.Get("GRMHD")->Param<std::string>("fname_fill1");
    auto fname_fill2 = pmb->packages.Get("GRMHD")->Param<std::string>("fname_fill2");
    auto fname_fill3 = pmb->packages.Get("GRMHD")->Param<std::string>("fname_fill3");
    auto fname_fill4 = pmb->packages.Get("GRMHD")->Param<std::string>("fname_fill4");
    auto fname_fill5 = pmb->packages.Get("GRMHD")->Param<std::string>("fname_fill5");
    auto fname_fill6 = pmb->packages.Get("GRMHD")->Param<std::string>("fname_fill6");
    auto fname_fill7 = pmb->packages.Get("GRMHD")->Param<std::string>("fname_fill7");
    const bool should_fill1 = !(fname_fill1 == "none");
    const bool should_fill2 = !(fname_fill2 == "none");
    const bool should_fill3 = !(fname_fill3 == "none");
    const bool should_fill4 = !(fname_fill4 == "none");
    const bool should_fill5 = !(fname_fill5 == "none");
    const bool should_fill6 = !(fname_fill6 == "none");
    const bool should_fill7 = !(fname_fill7 == "none");
    const Real fx1min = pmb->packages.Get("GRMHD")->Param<Real>("rx1min");
    const Real fx1max = pmb->packages.Get("GRMHD")->Param<Real>("rx1max");
    const Real fx1min_f1 = pmb->packages.Get("GRMHD")->Param<Real>("rx1min_f1");
    const Real fx1max_f1 = pmb->packages.Get("GRMHD")->Param<Real>("rx1max_f1");
    const Real fx1min_f2 = pmb->packages.Get("GRMHD")->Param<Real>("rx1min_f2");
    const Real fx1max_f2 = pmb->packages.Get("GRMHD")->Param<Real>("rx1max_f2");
    const Real fx1min_f3 = pmb->packages.Get("GRMHD")->Param<Real>("rx1min_f3");
    const Real fx1max_f3 = pmb->packages.Get("GRMHD")->Param<Real>("rx1max_f3");
    const Real fx1min_f4 = pmb->packages.Get("GRMHD")->Param<Real>("rx1min_f4");
    const Real fx1max_f4 = pmb->packages.Get("GRMHD")->Param<Real>("rx1max_f4");
    const Real fx1min_f5 = pmb->packages.Get("GRMHD")->Param<Real>("rx1min_f5");
    const Real fx1max_f5 = pmb->packages.Get("GRMHD")->Param<Real>("rx1max_f5");
    const Real fx1min_f6 = pmb->packages.Get("GRMHD")->Param<Real>("rx1min_f6");
    const Real fx1max_f6 = pmb->packages.Get("GRMHD")->Param<Real>("rx1max_f6");
    const Real fx1min_f7 = pmb->packages.Get("GRMHD")->Param<Real>("rx1min_f7");
    const Real fx1max_f7 = pmb->packages.Get("GRMHD")->Param<Real>("rx1max_f7");
    //const Real dx1 = (fx1max - fx1min) / n1tot;
    const Real x2min = pmb->packages.Get("GRMHD")->Param<Real>("x2min");
    const Real x2max = pmb->packages.Get("GRMHD")->Param<Real>("x2max");
    const Real x3min = pmb->packages.Get("GRMHD")->Param<Real>("x3min");
    const Real x3max = pmb->packages.Get("GRMHD")->Param<Real>("x3max");
    const GReal dx[GR_DIM] = {0., (fx1max - fx1min)/n1tot,
                                  (x2max - x2min)/n2tot,
                                  (x3max - x3min)/n3tot};
    const bool fghostzones = pmb->packages.Get("GRMHD")->Param<bool>("rghostzones");
    int fnghost = pmb->packages.Get("GRMHD")->Param<int>("rnghost");
    const Real fx1min_ghost = fx1min - fnghost*dx[1];
    const Real fx1max_ghost = fx1max + fnghost*dx[1];
    PackIndexMap prims_map, cons_map;
    auto P = GRMHD::PackMHDPrims(rc, prims_map);
    auto U = GRMHD::PackMHDCons(rc, cons_map);
    const VarMap m_u(cons_map, true), m_p(prims_map, false);
    
    if ((domain != IndexDomain::outer_x1) && (domain != IndexDomain::inner_x1)) { 
        // read from a restart file and save it to static GridScalar
        //cout << "Hyerin: reading files" << endl;


        if (! fghostzones) fnghost=0; // reset to 0
        int x3factor=1;
        if (n3tot <= 1) x3factor=0; // if less than 3D, do not add ghosts in x3
        hsize_t length[GR_DIM] = {nBlocks,
                                    n1mb+2*fnghost,
                                    n2mb+2*fnghost,
                                    n3mb+2*fnghost*x3factor}; 
        const int block_sz = length[0]*length[1]*length[2]*length[3];
        hsize_t f_length[GR_DIM] = {f_nBlocks,
                                    f_n1mb+2*fnghost,
                                    f_n2mb+2*fnghost,
                                    f_n3mb+2*fnghost*x3factor}; 
        const int f_block_sz = f_length[0]*f_length[1]*f_length[2]*f_length[3];
        //std::cout << "lengths " << length[0]  << " " << length[1] <<" " <<  length[2]<<" " << length[3] << std::endl;
        //printf("lengths %i %i %i %i \n", length[0], length[1], length[2], length[3]);
        
        
        // read from file and stored in device Hyerin (10/18/2022)
        GridScalar x1_f_device("x1_f_device", length[0], length[1]); 
        GridScalar x2_f_device("x2_f_device", length[0], length[2]); 
        GridScalar x3_f_device("x3_f_device", length[0], length[3]); 
        GridScalar rho_f_device("rho_f_device", length[0], length[3], length[2], length[1]); 
        GridScalar u_f_device("u_f_device", length[0], length[3], length[2], length[1]); 
        GridVector uvec_f_device("uvec_f_device", NVEC, length[0], length[3], length[2], length[1]); 
        GridVector B_f_device("B_f_device", NVEC, length[0], length[3], length[2], length[1]);
        auto x1_f_host = x1_f_device.GetHostMirror();
        auto x2_f_host = x2_f_device.GetHostMirror();
        auto x3_f_host = x3_f_device.GetHostMirror();
        auto rho_f_host = rho_f_device.GetHostMirror();
        auto u_f_host = u_f_device.GetHostMirror();
        auto uvec_f_host = uvec_f_device.GetHostMirror();
        auto B_f_host = B_f_device.GetHostMirror();
        // Hyerin (09/19/2022) : new attempt to read the file 
        hdf5_open(fname.c_str());
        hdf5_set_directory("/");
        Real *rho_file = new double[block_sz];
        Real *u_file = new double[block_sz];
        Real *uvec_file = new double[block_sz*3];
        Real *B_file = new double[block_sz*3];
        Real *x1_file = new double[length[0]*length[1]];
        Real *x2_file = new double[length[0]*length[2]];
        Real *x3_file = new double[length[0]*length[3]];
        //static hsize_t fdims[] = {length[0], length[3], length[2], length[1],1}; //outdated
        static hsize_t fdims[] = {length[0], 1, length[3], length[2], length[1]};
        //static hsize_t fdims_vec[] = {length[0], length[3], length[2], length[1],3}; //outdated
        static hsize_t fdims_vec[] = {length[0], 3, length[3], length[2], length[1]};
        static hsize_t fdims_x1[] = {length[0], length[1]};
        static hsize_t fdims_x2[] = {length[0], length[2]};
        static hsize_t fdims_x3[] = {length[0], length[3]};
        hsize_t fstart[] = {0, 0, 0, 0, 0};
        hsize_t fstart_x[] = {0, 0};
        hdf5_read_array(rho_file, "prims.rho", 5, fdims, fstart,fdims,fdims,fstart,H5T_IEEE_F64LE);
        hdf5_read_array(u_file, "prims.u", 5, fdims, fstart,fdims,fdims,fstart,H5T_IEEE_F64LE);
        hdf5_read_array(uvec_file, "prims.uvec", 5, fdims_vec, fstart,fdims_vec,fdims_vec,fstart,H5T_IEEE_F64LE);
        //if (include_B) hdf5_read_array(B_file, "prims.B", 5, fdims_vec, fstart,fdims_vec,fdims_vec,fstart,H5T_IEEE_F64LE);
        if (include_B) hdf5_read_array(B_file, "cons.B", 5, fdims_vec, fstart,fdims_vec,fdims_vec,fstart,H5T_IEEE_F64LE);
        hdf5_read_array(x1_file, "VolumeLocations/x", 2, fdims_x1, fstart_x,fdims_x1,fdims_x1,fstart_x,H5T_IEEE_F64LE);
        hdf5_read_array(x2_file, "VolumeLocations/y", 2, fdims_x2, fstart_x,fdims_x2,fdims_x2,fstart_x,H5T_IEEE_F64LE);
        hdf5_read_array(x3_file, "VolumeLocations/z", 2, fdims_x3, fstart_x,fdims_x3,fdims_x3,fstart_x,H5T_IEEE_F64LE);
        hdf5_close();
        
        GridScalar x1_fill1_device("x1_fill1_device", f_length[0], f_length[1]); 
        GridScalar x2_fill1_device("x2_fill1_device", f_length[0], f_length[2]); 
        GridScalar x3_fill1_device("x2_fill1_device", f_length[0], f_length[3]); 
        GridScalar rho_fill1_device("rho_fill1_device", f_length[0], f_length[3], f_length[2], f_length[1]); 
        GridScalar u_fill1_device("u_fill1_device", f_length[0], f_length[3], f_length[2], f_length[1]); 
        GridVector uvec_fill1_device("uvec_fill1_device", NVEC, f_length[0], f_length[3], f_length[2], f_length[1]); 
        GridVector B_fill1_device("B_fill1_device", NVEC, f_length[0], f_length[3], f_length[2], f_length[1]); 
        auto x1_fill1_host = x1_fill1_device.GetHostMirror();
        auto x2_fill1_host = x2_fill1_device.GetHostMirror();
        auto x3_fill1_host = x3_fill1_device.GetHostMirror();
        auto rho_fill1_host = rho_fill1_device.GetHostMirror();
        auto u_fill1_host = u_fill1_device.GetHostMirror();
        auto uvec_fill1_host = uvec_fill1_device.GetHostMirror();
        auto B_fill1_host = B_fill1_device.GetHostMirror();
        Real *rho_filefill1 = new double[f_block_sz];
        Real *u_filefill1 = new double[f_block_sz];
        Real *uvec_filefill1 = new double[f_block_sz*3];
        Real *B_filefill1 = new double[f_block_sz*3];
        Real *x1_filefill1 = new double[f_length[0]*f_length[1]];
        Real *x2_filefill1 = new double[f_length[0]*f_length[2]];
        Real *x3_filefill1 = new double[f_length[0]*f_length[3]];
        static hsize_t fill_fdims[] = {f_length[0], 1, f_length[3], f_length[2], f_length[1]};
        static hsize_t fill_fdims_vec[] = {f_length[0], 3, f_length[3], f_length[2], f_length[1]};
        static hsize_t fill_fdims_x1[] = {f_length[0], f_length[1]};
        static hsize_t fill_fdims_x2[] = {f_length[0], f_length[2]};
        static hsize_t fill_fdims_x3[] = {f_length[0], f_length[3]};
        if (should_fill1) { // TODO: here I'm assuming fname and fname_fill has same dimensions, which is not always the case.
            hdf5_open(fname_fill1.c_str());
            hdf5_set_directory("/");
            hdf5_read_array(rho_filefill1, "prims.rho", 5, fill_fdims, fstart,fill_fdims,fill_fdims,fstart,H5T_IEEE_F64LE);
            hdf5_read_array(u_filefill1, "prims.u", 5, fill_fdims, fstart,fill_fdims,fill_fdims,fstart,H5T_IEEE_F64LE);
            hdf5_read_array(uvec_filefill1, "prims.uvec", 5, fill_fdims_vec, fstart,fill_fdims_vec,fill_fdims_vec,fstart,H5T_IEEE_F64LE);
            if (include_B) hdf5_read_array(B_filefill1, "cons.B", 5, fill_fdims_vec, fstart,fill_fdims_vec,fill_fdims_vec,fstart,H5T_IEEE_F64LE);
            hdf5_read_array(x1_filefill1, "VolumeLocations/x", 2, fill_fdims_x1, fstart_x,fill_fdims_x1,fill_fdims_x1,fstart_x,H5T_IEEE_F64LE);
            hdf5_read_array(x2_filefill1, "VolumeLocations/y", 2, fill_fdims_x2, fstart_x,fill_fdims_x2,fill_fdims_x2,fstart_x,H5T_IEEE_F64LE);
            hdf5_read_array(x3_filefill1, "VolumeLocations/z", 2, fill_fdims_x3, fstart_x,fill_fdims_x3,fill_fdims_x3,fstart_x,H5T_IEEE_F64LE);
            hdf5_close();
        }

        GridScalar x1_fill2_device("x1_fill2_device", f_length[0], f_length[1]); 
        GridScalar x2_fill2_device("x2_fill2_device", f_length[0], f_length[2]); 
        GridScalar x3_fill2_device("x2_fill2_device", f_length[0], f_length[3]); 
        GridScalar rho_fill2_device("rho_fill2_device", f_length[0], f_length[3], f_length[2], f_length[1]); 
        GridScalar u_fill2_device("u_fill2_device", f_length[0], f_length[3], f_length[2], f_length[1]); 
        GridVector uvec_fill2_device("uvec_fill2_device", NVEC, f_length[0], f_length[3], f_length[2], f_length[1]); 
        GridVector B_fill2_device("B_fill2_device", NVEC, f_length[0], f_length[3], f_length[2], f_length[1]); 
        auto x1_fill2_host = x1_fill2_device.GetHostMirror();
        auto x2_fill2_host = x2_fill2_device.GetHostMirror();
        auto x3_fill2_host = x3_fill2_device.GetHostMirror();
        auto rho_fill2_host = rho_fill2_device.GetHostMirror();
        auto u_fill2_host = u_fill2_device.GetHostMirror();
        auto uvec_fill2_host = uvec_fill2_device.GetHostMirror();
        auto B_fill2_host = B_fill2_device.GetHostMirror();
        Real *rho_filefill2 = new double[f_block_sz];
        Real *u_filefill2 = new double[f_block_sz];
        Real *uvec_filefill2 = new double[f_block_sz*3];
        Real *B_filefill2 = new double[f_block_sz*3];
        Real *x1_filefill2 = new double[f_length[0]*f_length[1]];
        Real *x2_filefill2 = new double[f_length[0]*f_length[2]];
        Real *x3_filefill2 = new double[f_length[0]*f_length[3]];
        if (should_fill2) { // TODO: here I'm assuming fname and fname_fill has same dimensions, which is not always the case.
            hdf5_open(fname_fill2.c_str());
            hdf5_set_directory("/");
            hdf5_read_array(rho_filefill2, "prims.rho", 5, fill_fdims, fstart,fill_fdims,fill_fdims,fstart,H5T_IEEE_F64LE);
            hdf5_read_array(u_filefill2, "prims.u", 5, fill_fdims, fstart,fill_fdims,fill_fdims,fstart,H5T_IEEE_F64LE);
            hdf5_read_array(uvec_filefill2, "prims.uvec", 5, fill_fdims_vec, fstart,fill_fdims_vec,fill_fdims_vec,fstart,H5T_IEEE_F64LE);
            if (include_B) hdf5_read_array(B_filefill2, "cons.B", 5, fill_fdims_vec, fstart,fill_fdims_vec,fill_fdims_vec,fstart,H5T_IEEE_F64LE);
            hdf5_read_array(x1_filefill2, "VolumeLocations/x", 2, fill_fdims_x1, fstart_x,fill_fdims_x1,fill_fdims_x1,fstart_x,H5T_IEEE_F64LE);
            hdf5_read_array(x2_filefill2, "VolumeLocations/y", 2, fill_fdims_x2, fstart_x,fill_fdims_x2,fill_fdims_x2,fstart_x,H5T_IEEE_F64LE);
            hdf5_read_array(x3_filefill2, "VolumeLocations/z", 2, fill_fdims_x3, fstart_x,fill_fdims_x3,fill_fdims_x3,fstart_x,H5T_IEEE_F64LE);
            hdf5_close();
        }

        GridScalar x1_fill3_device("x1_fill3_device", f_length[0], f_length[1]); 
        GridScalar x2_fill3_device("x2_fill3_device", f_length[0], f_length[2]); 
        GridScalar x3_fill3_device("x2_fill3_device", f_length[0], f_length[3]); 
        GridScalar rho_fill3_device("rho_fill3_device", f_length[0], f_length[3], f_length[2], f_length[1]); 
        GridScalar u_fill3_device("u_fill3_device", f_length[0], f_length[3], f_length[2], f_length[1]); 
        GridVector uvec_fill3_device("uvec_fill3_device", NVEC, f_length[0], f_length[3], f_length[2], f_length[1]); 
        GridVector B_fill3_device("B_fill3_device", NVEC, f_length[0], f_length[3], f_length[2], f_length[1]); 
        auto x1_fill3_host = x1_fill3_device.GetHostMirror();
        auto x2_fill3_host = x2_fill3_device.GetHostMirror();
        auto x3_fill3_host = x3_fill3_device.GetHostMirror();
        auto rho_fill3_host = rho_fill3_device.GetHostMirror();
        auto u_fill3_host = u_fill3_device.GetHostMirror();
        auto uvec_fill3_host = uvec_fill3_device.GetHostMirror();
        auto B_fill3_host = B_fill3_device.GetHostMirror();
        Real *rho_filefill3 = new double[f_block_sz];
        Real *u_filefill3 = new double[f_block_sz];
        Real *uvec_filefill3 = new double[f_block_sz*3];
        Real *B_filefill3 = new double[f_block_sz*3];
        Real *x1_filefill3 = new double[f_length[0]*f_length[1]];
        Real *x2_filefill3 = new double[f_length[0]*f_length[2]];
        Real *x3_filefill3 = new double[f_length[0]*f_length[3]];
        if (should_fill3) { // TODO: here I'm assuming fname and fname_fill has same dimensions, which is not always the case.
            hdf5_open(fname_fill3.c_str());
            hdf5_set_directory("/");
            hdf5_read_array(rho_filefill3, "prims.rho", 5, fill_fdims, fstart,fill_fdims,fill_fdims,fstart,H5T_IEEE_F64LE);
            hdf5_read_array(u_filefill3, "prims.u", 5, fill_fdims, fstart,fill_fdims,fill_fdims,fstart,H5T_IEEE_F64LE);
            hdf5_read_array(uvec_filefill3, "prims.uvec", 5, fill_fdims_vec, fstart,fill_fdims_vec,fill_fdims_vec,fstart,H5T_IEEE_F64LE);
            if (include_B) hdf5_read_array(B_filefill3, "cons.B", 5, fill_fdims_vec, fstart,fill_fdims_vec,fill_fdims_vec,fstart,H5T_IEEE_F64LE);
            hdf5_read_array(x1_filefill3, "VolumeLocations/x", 2, fill_fdims_x1, fstart_x,fill_fdims_x1,fill_fdims_x1,fstart_x,H5T_IEEE_F64LE);
            hdf5_read_array(x2_filefill3, "VolumeLocations/y", 2, fill_fdims_x2, fstart_x,fill_fdims_x2,fill_fdims_x2,fstart_x,H5T_IEEE_F64LE);
            hdf5_read_array(x3_filefill3, "VolumeLocations/z", 2, fill_fdims_x3, fstart_x,fill_fdims_x3,fill_fdims_x3,fstart_x,H5T_IEEE_F64LE);
            hdf5_close();
        }

        GridScalar x1_fill4_device("x1_fill4_device", f_length[0], f_length[1]); 
        GridScalar x2_fill4_device("x2_fill4_device", f_length[0], f_length[2]); 
        GridScalar x3_fill4_device("x2_fill4_device", f_length[0], f_length[3]); 
        GridScalar rho_fill4_device("rho_fill4_device", f_length[0], f_length[3], f_length[2], f_length[1]); 
        GridScalar u_fill4_device("u_fill4_device", f_length[0], f_length[3], f_length[2], f_length[1]); 
        GridVector uvec_fill4_device("uvec_fill4_device", NVEC, f_length[0], f_length[3], f_length[2], f_length[1]); 
        GridVector B_fill4_device("B_fill4_device", NVEC, f_length[0], f_length[3], f_length[2], f_length[1]); 
        auto x1_fill4_host = x1_fill4_device.GetHostMirror();
        auto x2_fill4_host = x2_fill4_device.GetHostMirror();
        auto x3_fill4_host = x3_fill4_device.GetHostMirror();
        auto rho_fill4_host = rho_fill4_device.GetHostMirror();
        auto u_fill4_host = u_fill4_device.GetHostMirror();
        auto uvec_fill4_host = uvec_fill4_device.GetHostMirror();
        auto B_fill4_host = B_fill4_device.GetHostMirror();
        Real *rho_filefill4 = new double[f_block_sz];
        Real *u_filefill4 = new double[f_block_sz];
        Real *uvec_filefill4 = new double[f_block_sz*3];
        Real *B_filefill4 = new double[f_block_sz*3];
        Real *x1_filefill4 = new double[f_length[0]*f_length[1]];
        Real *x2_filefill4 = new double[f_length[0]*f_length[2]];
        Real *x3_filefill4 = new double[f_length[0]*f_length[3]];
        if (should_fill4) { // TODO: here I'm assuming fname and fname_fill has same dimensions, which is not always the case.
            hdf5_open(fname_fill4.c_str());
            hdf5_set_directory("/");
            hdf5_read_array(rho_filefill4, "prims.rho", 5, fill_fdims, fstart,fill_fdims,fill_fdims,fstart,H5T_IEEE_F64LE);
            hdf5_read_array(u_filefill4, "prims.u", 5, fill_fdims, fstart,fill_fdims,fill_fdims,fstart,H5T_IEEE_F64LE);
            hdf5_read_array(uvec_filefill4, "prims.uvec", 5, fill_fdims_vec, fstart,fill_fdims_vec,fill_fdims_vec,fstart,H5T_IEEE_F64LE);
            if (include_B) hdf5_read_array(B_filefill4, "cons.B", 5, fill_fdims_vec, fstart,fill_fdims_vec,fill_fdims_vec,fstart,H5T_IEEE_F64LE);
            hdf5_read_array(x1_filefill4, "VolumeLocations/x", 2, fill_fdims_x1, fstart_x,fill_fdims_x1,fill_fdims_x1,fstart_x,H5T_IEEE_F64LE);
            hdf5_read_array(x2_filefill4, "VolumeLocations/y", 2, fill_fdims_x2, fstart_x,fill_fdims_x2,fill_fdims_x2,fstart_x,H5T_IEEE_F64LE);
            hdf5_read_array(x3_filefill4, "VolumeLocations/z", 2, fill_fdims_x3, fstart_x,fill_fdims_x3,fill_fdims_x3,fstart_x,H5T_IEEE_F64LE);
            hdf5_close();
        }

        GridScalar x1_fill5_device("x1_fill5_device", f_length[0], f_length[1]); 
        GridScalar x2_fill5_device("x2_fill5_device", f_length[0], f_length[2]); 
        GridScalar x3_fill5_device("x2_fill5_device", f_length[0], f_length[3]); 
        GridScalar rho_fill5_device("rho_fill5_device", f_length[0], f_length[3], f_length[2], f_length[1]); 
        GridScalar u_fill5_device("u_fill5_device", f_length[0], f_length[3], f_length[2], f_length[1]); 
        GridVector uvec_fill5_device("uvec_fill5_device", NVEC, f_length[0], f_length[3], f_length[2], f_length[1]); 
        GridVector B_fill5_device("B_fill5_device", NVEC, f_length[0], f_length[3], f_length[2], f_length[1]); 
        auto x1_fill5_host = x1_fill5_device.GetHostMirror();
        auto x2_fill5_host = x2_fill5_device.GetHostMirror();
        auto x3_fill5_host = x3_fill5_device.GetHostMirror();
        auto rho_fill5_host = rho_fill5_device.GetHostMirror();
        auto u_fill5_host = u_fill5_device.GetHostMirror();
        auto uvec_fill5_host = uvec_fill5_device.GetHostMirror();
        auto B_fill5_host = B_fill5_device.GetHostMirror();
        Real *rho_filefill5 = new double[f_block_sz];
        Real *u_filefill5 = new double[f_block_sz];
        Real *uvec_filefill5 = new double[f_block_sz*3];
        Real *B_filefill5 = new double[f_block_sz*3];
        Real *x1_filefill5 = new double[f_length[0]*f_length[1]];
        Real *x2_filefill5 = new double[f_length[0]*f_length[2]];
        Real *x3_filefill5 = new double[f_length[0]*f_length[3]];
        if (should_fill5) { // TODO: here I'm assuming fname and fname_fill has same dimensions, which is not always the case.
            hdf5_open(fname_fill5.c_str());
            hdf5_set_directory("/");
            hdf5_read_array(rho_filefill5, "prims.rho", 5, fill_fdims, fstart,fill_fdims,fill_fdims,fstart,H5T_IEEE_F64LE);
            hdf5_read_array(u_filefill5, "prims.u", 5, fill_fdims, fstart,fill_fdims,fill_fdims,fstart,H5T_IEEE_F64LE);
            hdf5_read_array(uvec_filefill5, "prims.uvec", 5, fill_fdims_vec, fstart,fill_fdims_vec,fill_fdims_vec,fstart,H5T_IEEE_F64LE);
            if (include_B) hdf5_read_array(B_filefill5, "cons.B", 5, fill_fdims_vec, fstart,fill_fdims_vec,fill_fdims_vec,fstart,H5T_IEEE_F64LE);
            hdf5_read_array(x1_filefill5, "VolumeLocations/x", 2, fill_fdims_x1, fstart_x,fill_fdims_x1,fill_fdims_x1,fstart_x,H5T_IEEE_F64LE);
            hdf5_read_array(x2_filefill5, "VolumeLocations/y", 2, fill_fdims_x2, fstart_x,fill_fdims_x2,fill_fdims_x2,fstart_x,H5T_IEEE_F64LE);
            hdf5_read_array(x3_filefill5, "VolumeLocations/z", 2, fill_fdims_x3, fstart_x,fill_fdims_x3,fill_fdims_x3,fstart_x,H5T_IEEE_F64LE);
            hdf5_close();
        }

        GridScalar x1_fill6_device("x1_fill6_device", f_length[0], f_length[1]); 
        GridScalar x2_fill6_device("x2_fill6_device", f_length[0], f_length[2]); 
        GridScalar x3_fill6_device("x2_fill6_device", f_length[0], f_length[3]); 
        GridScalar rho_fill6_device("rho_fill6_device", f_length[0], f_length[3], f_length[2], f_length[1]); 
        GridScalar u_fill6_device("u_fill6_device", f_length[0], f_length[3], f_length[2], f_length[1]); 
        GridVector uvec_fill6_device("uvec_fill6_device", NVEC, f_length[0], f_length[3], f_length[2], f_length[1]); 
        GridVector B_fill6_device("B_fill6_device", NVEC, f_length[0], f_length[3], f_length[2], f_length[1]); 
        auto x1_fill6_host = x1_fill6_device.GetHostMirror();
        auto x2_fill6_host = x2_fill6_device.GetHostMirror();
        auto x3_fill6_host = x3_fill6_device.GetHostMirror();
        auto rho_fill6_host = rho_fill6_device.GetHostMirror();
        auto u_fill6_host = u_fill6_device.GetHostMirror();
        auto uvec_fill6_host = uvec_fill6_device.GetHostMirror();
        auto B_fill6_host = B_fill6_device.GetHostMirror();
        Real *rho_filefill6 = new double[f_block_sz];
        Real *u_filefill6 = new double[f_block_sz];
        Real *uvec_filefill6 = new double[f_block_sz*3];
        Real *B_filefill6 = new double[f_block_sz*3];
        Real *x1_filefill6 = new double[f_length[0]*f_length[1]];
        Real *x2_filefill6 = new double[f_length[0]*f_length[2]];
        Real *x3_filefill6 = new double[f_length[0]*f_length[3]];
        if (should_fill6) { // TODO: here I'm assuming fname and fname_fill has same dimensions, which is not always the case.
            hdf5_open(fname_fill6.c_str());
            hdf5_set_directory("/");
            hdf5_read_array(rho_filefill6, "prims.rho", 5, fill_fdims, fstart,fill_fdims,fill_fdims,fstart,H5T_IEEE_F64LE);
            hdf5_read_array(u_filefill6, "prims.u", 5, fill_fdims, fstart,fill_fdims,fill_fdims,fstart,H5T_IEEE_F64LE);
            hdf5_read_array(uvec_filefill6, "prims.uvec", 5, fill_fdims_vec, fstart,fill_fdims_vec,fill_fdims_vec,fstart,H5T_IEEE_F64LE);
            if (include_B) hdf5_read_array(B_filefill6, "cons.B", 5, fill_fdims_vec, fstart,fill_fdims_vec,fill_fdims_vec,fstart,H5T_IEEE_F64LE);
            hdf5_read_array(x1_filefill6, "VolumeLocations/x", 2, fill_fdims_x1, fstart_x,fill_fdims_x1,fill_fdims_x1,fstart_x,H5T_IEEE_F64LE);
            hdf5_read_array(x2_filefill6, "VolumeLocations/y", 2, fill_fdims_x2, fstart_x,fill_fdims_x2,fill_fdims_x2,fstart_x,H5T_IEEE_F64LE);
            hdf5_read_array(x3_filefill6, "VolumeLocations/z", 2, fill_fdims_x3, fstart_x,fill_fdims_x3,fill_fdims_x3,fstart_x,H5T_IEEE_F64LE);
            hdf5_close();
        }

        GridScalar x1_fill7_device("x1_fill7_device", f_length[0], f_length[1]); 
        GridScalar x2_fill7_device("x2_fill7_device", f_length[0], f_length[2]); 
        GridScalar x3_fill7_device("x2_fill7_device", f_length[0], f_length[3]); 
        GridScalar rho_fill7_device("rho_fill7_device", f_length[0], f_length[3], f_length[2], f_length[1]); 
        GridScalar u_fill7_device("u_fill7_device", f_length[0], f_length[3], f_length[2], f_length[1]); 
        GridVector uvec_fill7_device("uvec_fill7_device", NVEC, f_length[0], f_length[3], f_length[2], f_length[1]); 
        GridVector B_fill7_device("B_fill7_device", NVEC, f_length[0], f_length[3], f_length[2], f_length[1]); 
        auto x1_fill7_host = x1_fill7_device.GetHostMirror();
        auto x2_fill7_host = x2_fill7_device.GetHostMirror();
        auto x3_fill7_host = x3_fill7_device.GetHostMirror();
        auto rho_fill7_host = rho_fill7_device.GetHostMirror();
        auto u_fill7_host = u_fill7_device.GetHostMirror();
        auto uvec_fill7_host = uvec_fill7_device.GetHostMirror();
        auto B_fill7_host = B_fill7_device.GetHostMirror();
        Real *rho_filefill7 = new double[f_block_sz];
        Real *u_filefill7 = new double[f_block_sz];
        Real *uvec_filefill7 = new double[f_block_sz*3];
        Real *B_filefill7 = new double[f_block_sz*3];
        Real *x1_filefill7 = new double[f_length[0]*f_length[1]];
        Real *x2_filefill7 = new double[f_length[0]*f_length[2]];
        Real *x3_filefill7 = new double[f_length[0]*f_length[3]];
        if (should_fill7) { // TODO: here I'm assuming fname and fname_fill has same dimensions, which is not always the case.
            hdf5_open(fname_fill7.c_str());
            hdf5_set_directory("/");
            hdf5_read_array(rho_filefill7, "prims.rho", 5, fill_fdims, fstart,fill_fdims,fill_fdims,fstart,H5T_IEEE_F64LE);
            hdf5_read_array(u_filefill7, "prims.u", 5, fill_fdims, fstart,fill_fdims,fill_fdims,fstart,H5T_IEEE_F64LE);
            hdf5_read_array(uvec_filefill7, "prims.uvec", 5, fill_fdims_vec, fstart,fill_fdims_vec,fill_fdims_vec,fstart,H5T_IEEE_F64LE);
            if (include_B) hdf5_read_array(B_filefill7, "cons.B", 5, fill_fdims_vec, fstart,fill_fdims_vec,fill_fdims_vec,fstart,H5T_IEEE_F64LE);
            hdf5_read_array(x1_filefill7, "VolumeLocations/x", 2, fill_fdims_x1, fstart_x,fill_fdims_x1,fill_fdims_x1,fstart_x,H5T_IEEE_F64LE);
            hdf5_read_array(x2_filefill7, "VolumeLocations/y", 2, fill_fdims_x2, fstart_x,fill_fdims_x2,fill_fdims_x2,fstart_x,H5T_IEEE_F64LE);
            hdf5_read_array(x3_filefill7, "VolumeLocations/z", 2, fill_fdims_x3, fstart_x,fill_fdims_x3,fill_fdims_x3,fstart_x,H5T_IEEE_F64LE);
            hdf5_close();
        }

        // save the grid coordinate values to host array
        for (int iblocktemp = 0; iblocktemp < length[0]; iblocktemp++) {
            for (int itemp = 0; itemp < length[1]; itemp++) {
                x1_f_host(iblocktemp,itemp) = x1_file[length[1]*iblocktemp+itemp];
            } for (int jtemp = 0; jtemp < length[2]; jtemp++) {
                x2_f_host(iblocktemp,jtemp) = x2_file[length[2]*iblocktemp+jtemp];
            } for (int ktemp = 0; ktemp < length[3]; ktemp++) {
                x3_f_host(iblocktemp,ktemp) = x3_file[length[3]*iblocktemp+ktemp];
            }
        }
        for (int iblocktemp = 0; iblocktemp < f_length[0]; iblocktemp++) {
            for (int itemp = 0; itemp < f_length[1]; itemp++) {
                if (should_fill1) x1_fill1_host(iblocktemp,itemp) = x1_filefill1[f_length[1]*iblocktemp+itemp];
                if (should_fill2) x1_fill2_host(iblocktemp,itemp) = x1_filefill2[f_length[1]*iblocktemp+itemp];
                if (should_fill3) x1_fill3_host(iblocktemp,itemp) = x1_filefill3[f_length[1]*iblocktemp+itemp];
                if (should_fill4) x1_fill4_host(iblocktemp,itemp) = x1_filefill4[f_length[1]*iblocktemp+itemp];
                if (should_fill5) x1_fill5_host(iblocktemp,itemp) = x1_filefill5[f_length[1]*iblocktemp+itemp];
                if (should_fill6) x1_fill6_host(iblocktemp,itemp) = x1_filefill6[f_length[1]*iblocktemp+itemp];
                if (should_fill7) x1_fill7_host(iblocktemp,itemp) = x1_filefill7[f_length[1]*iblocktemp+itemp];
            } for (int jtemp = 0; jtemp < f_length[2]; jtemp++) {
                if (should_fill1) x2_fill1_host(iblocktemp,jtemp) = x2_filefill1[f_length[2]*iblocktemp+jtemp];
                if (should_fill2) x2_fill2_host(iblocktemp,jtemp) = x2_filefill2[f_length[2]*iblocktemp+jtemp];
                if (should_fill3) x2_fill3_host(iblocktemp,jtemp) = x2_filefill3[f_length[2]*iblocktemp+jtemp];
                if (should_fill4) x2_fill4_host(iblocktemp,jtemp) = x2_filefill4[f_length[2]*iblocktemp+jtemp];
                if (should_fill5) x2_fill5_host(iblocktemp,jtemp) = x2_filefill5[f_length[2]*iblocktemp+jtemp];
                if (should_fill6) x2_fill6_host(iblocktemp,jtemp) = x2_filefill6[f_length[2]*iblocktemp+jtemp];
                if (should_fill7) x2_fill7_host(iblocktemp,jtemp) = x2_filefill7[f_length[2]*iblocktemp+jtemp];
            } for (int ktemp = 0; ktemp < f_length[3]; ktemp++) {
                if (should_fill1) x3_fill1_host(iblocktemp,ktemp) = x3_filefill1[f_length[3]*iblocktemp+ktemp];
                if (should_fill2) x3_fill2_host(iblocktemp,ktemp) = x3_filefill2[f_length[3]*iblocktemp+ktemp];
                if (should_fill3) x3_fill3_host(iblocktemp,ktemp) = x3_filefill3[f_length[3]*iblocktemp+ktemp];
                if (should_fill4) x3_fill4_host(iblocktemp,ktemp) = x3_filefill4[f_length[3]*iblocktemp+ktemp];
                if (should_fill5) x3_fill5_host(iblocktemp,ktemp) = x3_filefill5[f_length[3]*iblocktemp+ktemp];
                if (should_fill6) x3_fill6_host(iblocktemp,ktemp) = x3_filefill6[f_length[3]*iblocktemp+ktemp];
                if (should_fill7) x3_fill7_host(iblocktemp,ktemp) = x3_filefill7[f_length[3]*iblocktemp+ktemp];
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
                            vector_file_index = length[1]*(length[2]*(length[3]*(3*iblocktemp+ltemp)+ktemp)+jtemp)+itemp;
                            
                            uvec_f_host(ltemp,iblocktemp,ktemp,jtemp,itemp) = uvec_file[vector_file_index];
                            if (include_B) B_f_host(ltemp,iblocktemp,ktemp,jtemp,itemp) = B_file[vector_file_index];
                        }
                    }
                }
            }
        }
        for (int iblocktemp = 0; iblocktemp < f_length[0]; iblocktemp++) {
            for (int itemp = 0; itemp < f_length[1]; itemp++) {
                for (int jtemp = 0; jtemp < f_length[2]; jtemp++) {
                    for (int ktemp = 0; ktemp < f_length[3]; ktemp++) {
                        scalar_file_index = f_length[1]*(f_length[2]*(f_length[3]*iblocktemp+ktemp)+jtemp)+itemp;
                        if (should_fill1) {
                            rho_fill1_host(iblocktemp,ktemp,jtemp,itemp) = rho_filefill1[scalar_file_index];
                            u_fill1_host(iblocktemp,ktemp,jtemp,itemp) = u_filefill1[scalar_file_index];
                        }
                        if (should_fill2) {
                            rho_fill2_host(iblocktemp,ktemp,jtemp,itemp) = rho_filefill2[scalar_file_index];
                            u_fill2_host(iblocktemp,ktemp,jtemp,itemp) = u_filefill2[scalar_file_index];
                        }
                        if (should_fill3) {
                            rho_fill3_host(iblocktemp,ktemp,jtemp,itemp) = rho_filefill3[scalar_file_index];
                            u_fill3_host(iblocktemp,ktemp,jtemp,itemp) = u_filefill3[scalar_file_index];
                        }
                        if (should_fill4) {
                            rho_fill4_host(iblocktemp,ktemp,jtemp,itemp) = rho_filefill4[scalar_file_index];
                            u_fill4_host(iblocktemp,ktemp,jtemp,itemp) = u_filefill4[scalar_file_index];
                        }
                        if (should_fill5) {
                            rho_fill5_host(iblocktemp,ktemp,jtemp,itemp) = rho_filefill5[scalar_file_index];
                            u_fill5_host(iblocktemp,ktemp,jtemp,itemp) = u_filefill5[scalar_file_index];
                        }
                        if (should_fill6) {
                            rho_fill6_host(iblocktemp,ktemp,jtemp,itemp) = rho_filefill6[scalar_file_index];
                            u_fill6_host(iblocktemp,ktemp,jtemp,itemp) = u_filefill6[scalar_file_index];
                        }
                        if (should_fill7) {
                            rho_fill7_host(iblocktemp,ktemp,jtemp,itemp) = rho_filefill7[scalar_file_index];
                            u_fill7_host(iblocktemp,ktemp,jtemp,itemp) = u_filefill7[scalar_file_index];
                        }
                        for (int ltemp = 0; ltemp < 3; ltemp++) {
                            vector_file_index = f_length[1]*(f_length[2]*(f_length[3]*(3*iblocktemp+ltemp)+ktemp)+jtemp)+itemp;
                            if (should_fill1) {
                                uvec_fill1_host(ltemp,iblocktemp,ktemp,jtemp,itemp) = uvec_filefill1[vector_file_index];
                                if (include_B) B_fill1_host(ltemp,iblocktemp,ktemp,jtemp,itemp) = B_filefill1[vector_file_index];
                            }
                            if (should_fill2) {
                                uvec_fill2_host(ltemp,iblocktemp,ktemp,jtemp,itemp) = uvec_filefill2[vector_file_index];
                                if (include_B) B_fill2_host(ltemp,iblocktemp,ktemp,jtemp,itemp) = B_filefill2[vector_file_index];
                            }
                            if (should_fill3) {
                                uvec_fill3_host(ltemp,iblocktemp,ktemp,jtemp,itemp) = uvec_filefill3[vector_file_index];
                                if (include_B) B_fill3_host(ltemp,iblocktemp,ktemp,jtemp,itemp) = B_filefill3[vector_file_index];
                            }
                            if (should_fill4) {
                                uvec_fill4_host(ltemp,iblocktemp,ktemp,jtemp,itemp) = uvec_filefill4[vector_file_index];
                                if (include_B) B_fill4_host(ltemp,iblocktemp,ktemp,jtemp,itemp) = B_filefill4[vector_file_index];
                            }
                            if (should_fill5) {
                                uvec_fill5_host(ltemp,iblocktemp,ktemp,jtemp,itemp) = uvec_filefill5[vector_file_index];
                                if (include_B) B_fill5_host(ltemp,iblocktemp,ktemp,jtemp,itemp) = B_filefill5[vector_file_index];
                            }
                            if (should_fill6) {
                                uvec_fill6_host(ltemp,iblocktemp,ktemp,jtemp,itemp) = uvec_filefill6[vector_file_index];
                                if (include_B) B_fill6_host(ltemp,iblocktemp,ktemp,jtemp,itemp) = B_filefill6[vector_file_index];
                            }
                            if (should_fill7) {
                                uvec_fill7_host(ltemp,iblocktemp,ktemp,jtemp,itemp) = uvec_filefill7[vector_file_index];
                                if (include_B) B_fill7_host(ltemp,iblocktemp,ktemp,jtemp,itemp) = B_filefill7[vector_file_index];
                            }
                        }
                    }
                }
            }
        }
        //std::cout << "Hyerin: first five Bs" << B_file[0] << " " << B_file[1] << " " << B_file[2] << " " << B_file[3] << " " << B_file[4] << std::endl; 
        //std::cout << "Hyerin: 6,7,8,9,10 B_f " << B_f_host(0,0,0,0,6) << " " << B_f_host(0,0,0,0,7) << " " << B_f_host(0,0,0,0,8) << " " << B_f_host(0,0,0,0,9) << " " << B_f_host(0,0,0,0,10) << std::endl; 

      
        // Deep copy to device
        x1_f_device.DeepCopy(x1_f_host);
        x2_f_device.DeepCopy(x2_f_host);
        x3_f_device.DeepCopy(x3_f_host);
        rho_f_device.DeepCopy(rho_f_host);
        u_f_device.DeepCopy(u_f_host);
        uvec_f_device.DeepCopy(uvec_f_host);
        if (include_B) B_f_device.DeepCopy(B_f_host);
        if (should_fill1) {
            x1_fill1_device.DeepCopy(x1_fill1_host);
            x2_fill1_device.DeepCopy(x2_fill1_host);
            x3_fill1_device.DeepCopy(x3_fill1_host);
            rho_fill1_device.DeepCopy(rho_fill1_host);
            u_fill1_device.DeepCopy(u_fill1_host);
            uvec_fill1_device.DeepCopy(uvec_fill1_host);
            if (include_B) B_fill1_device.DeepCopy(B_fill1_host);
        }
        if (should_fill2) {
            x1_fill2_device.DeepCopy(x1_fill2_host);
            x2_fill2_device.DeepCopy(x2_fill2_host);
            x3_fill2_device.DeepCopy(x3_fill2_host);
            rho_fill2_device.DeepCopy(rho_fill2_host);
            u_fill2_device.DeepCopy(u_fill2_host);
            uvec_fill2_device.DeepCopy(uvec_fill2_host);
            if (include_B) B_fill2_device.DeepCopy(B_fill2_host);
        }
        if (should_fill3) {
            x1_fill3_device.DeepCopy(x1_fill3_host);
            x2_fill3_device.DeepCopy(x2_fill3_host);
            x3_fill3_device.DeepCopy(x3_fill3_host);
            rho_fill3_device.DeepCopy(rho_fill3_host);
            u_fill3_device.DeepCopy(u_fill3_host);
            uvec_fill3_device.DeepCopy(uvec_fill3_host);
            if (include_B) B_fill3_device.DeepCopy(B_fill3_host);
        }
        if (should_fill4) {
            x1_fill4_device.DeepCopy(x1_fill4_host);
            x2_fill4_device.DeepCopy(x2_fill4_host);
            x3_fill4_device.DeepCopy(x3_fill4_host);
            rho_fill4_device.DeepCopy(rho_fill4_host);
            u_fill4_device.DeepCopy(u_fill4_host);
            uvec_fill4_device.DeepCopy(uvec_fill4_host);
            if (include_B) B_fill4_device.DeepCopy(B_fill4_host);
        }
        if (should_fill5) {
            x1_fill5_device.DeepCopy(x1_fill5_host);
            x2_fill5_device.DeepCopy(x2_fill5_host);
            x3_fill5_device.DeepCopy(x3_fill5_host);
            rho_fill5_device.DeepCopy(rho_fill5_host);
            u_fill5_device.DeepCopy(u_fill5_host);
            uvec_fill5_device.DeepCopy(uvec_fill5_host);
            if (include_B) B_fill5_device.DeepCopy(B_fill5_host);
        }
        if (should_fill6) {
            x1_fill6_device.DeepCopy(x1_fill6_host);
            x2_fill6_device.DeepCopy(x2_fill6_host);
            x3_fill6_device.DeepCopy(x3_fill6_host);
            rho_fill6_device.DeepCopy(rho_fill6_host);
            u_fill6_device.DeepCopy(u_fill6_host);
            uvec_fill6_device.DeepCopy(uvec_fill6_host);
            if (include_B) B_fill6_device.DeepCopy(B_fill6_host);
        }
        if (should_fill7) {
            x1_fill7_device.DeepCopy(x1_fill7_host);
            x2_fill7_device.DeepCopy(x2_fill7_host);
            x3_fill7_device.DeepCopy(x3_fill7_host);
            rho_fill7_device.DeepCopy(rho_fill7_host);
            u_fill7_device.DeepCopy(u_fill7_host);
            uvec_fill7_device.DeepCopy(uvec_fill7_host);
            if (include_B) B_fill7_device.DeepCopy(B_fill7_host);
        }
        //if (pin->GetOrAddString("b_field", "type", "none") != "none") {
        //    B_P.DeepCopy(B_host);
        //}
        Kokkos::fence();

        // Host-side interpolate & copy into the mirror array
        pmb->par_for("copy_restart_state_kharma", ks, ke, js, je, is, ie,
            KOKKOS_LAMBDA_3D {
                GReal X[GR_DIM];
                G.coord(k, j, i, Loci::center, X);
                if ((should_fill1) && (X[1]>=fx1min_f1 - fnghost*dx[1]) && (X[1]<=fx1max_f1 + fnghost*dx[1])) { // fill with the fname_fill1
                    get_prim_restart_kharma(X, P, m_p, dx, f_length,
                        x1_fill1_device, x2_fill1_device, x3_fill1_device, rho_fill1_device, u_fill1_device, uvec_fill1_device,
                        k, j, i);
                    if (include_B) get_B_restart_kharma(X, P, m_p, dx, f_length, x1_fill1_device, x2_fill1_device, x3_fill1_device, B_fill1_device, B_Save, k, j, i);
                // TODO Hyerin (07/23/23) should also use x1min x1max including ghosts?
                } else if ((should_fill2) && (X[1]>=fx1min_f2 - fnghost*dx[1]) && (X[1]<=fx1max_f2 + fnghost*dx[1])) { // fill with the fname_fill2
                    get_prim_restart_kharma(X, P, m_p, dx, f_length,
                        x1_fill2_device, x2_fill2_device, x3_fill2_device, rho_fill2_device, u_fill2_device, uvec_fill2_device,
                        k, j, i);
                    if (include_B) get_B_restart_kharma(X, P, m_p, dx, f_length, x1_fill2_device, x2_fill2_device, x3_fill2_device, B_fill2_device, B_Save, k, j, i);
                } else if ((should_fill3) && (X[1]>=fx1min_f3) && (X[1]<=fx1max_f3)) { // fill with the fname_fill3
                    get_prim_restart_kharma(X, P, m_p, dx, f_length,
                        x1_fill3_device, x2_fill3_device, x3_fill3_device, rho_fill3_device, u_fill3_device, uvec_fill3_device,
                        k, j, i);
                    if (include_B) get_B_restart_kharma(X, P, m_p, dx, f_length, x1_fill3_device, x2_fill3_device, x3_fill3_device, B_fill3_device, B_Save, k, j, i);
                } else if ((should_fill4) && (X[1]>=fx1min_f4) && (X[1]<=fx1max_f4)) { // fill with the fname_fill4
                    get_prim_restart_kharma(X, P, m_p, dx, f_length,
                        x1_fill4_device, x2_fill4_device, x3_fill4_device, rho_fill4_device, u_fill4_device, uvec_fill4_device,
                        k, j, i);
                    if (include_B) get_B_restart_kharma(X, P, m_p, dx, f_length, x1_fill4_device, x2_fill4_device, x3_fill4_device, B_fill4_device, B_Save, k, j, i);
                } else if ((should_fill5) && (X[1]>=fx1min_f5) && (X[1]<=fx1max_f5)) { // fill with the fname_fill5
                    get_prim_restart_kharma(X, P, m_p, dx, f_length,
                        x1_fill5_device, x2_fill5_device, x3_fill5_device, rho_fill5_device, u_fill5_device, uvec_fill5_device,
                        k, j, i);
                    if (include_B) get_B_restart_kharma(X, P, m_p, dx, f_length, x1_fill5_device, x2_fill5_device, x3_fill5_device, B_fill5_device, B_Save, k, j, i);
                } else if ((should_fill6) && (X[1]>=fx1min_f6) && (X[1]<=fx1max_f6)) { // fill with the fname_fill6
                    get_prim_restart_kharma(X, P, m_p, dx, f_length,
                        x1_fill6_device, x2_fill6_device, x3_fill6_device, rho_fill6_device, u_fill6_device, uvec_fill6_device,
                        k, j, i);
                    if (include_B) get_B_restart_kharma(X, P, m_p, dx, f_length, x1_fill6_device, x2_fill6_device, x3_fill6_device, B_fill6_device, B_Save, k, j, i);
                } else if ((should_fill7) && (X[1]>=fx1min_f7) && (X[1]<=fx1max_f7)) { // fill with the fname_fill7
                    get_prim_restart_kharma(X, P, m_p, dx, f_length,
                        x1_fill7_device, x2_fill7_device, x3_fill7_device, rho_fill7_device, u_fill7_device, uvec_fill7_device,
                        k, j, i);
                    if (include_B) get_B_restart_kharma(X, P, m_p, dx, f_length, x1_fill7_device, x2_fill7_device, x3_fill7_device, B_fill7_device, B_Save, k, j, i);
                } else if ((X[1]>=fx1min_ghost) && (X[1]<=fx1max_ghost)) { // fill with the fname
                    get_prim_restart_kharma(X, P, m_p, dx, length,
                        x1_f_device, x2_f_device, x3_f_device, rho_f_device, u_f_device, uvec_f_device,
                        k, j, i);
                    if (include_B) get_B_restart_kharma(X, P, m_p, dx, length, x1_f_device, x2_f_device, x3_f_device, B_f_device, B_Save, k, j, i);
                } else {
                    printf("HYERIN: no corresponding file found to fill!\n");
                }
            }
        );
        //if (include_B) B_FluxCT::PtoU(rc,domain); // added for B fields
    }

   return TaskStatus::complete;
}
