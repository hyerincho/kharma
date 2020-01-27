/*
 * K/HARM -- Implementation of the HARM scheme for GRMHD made
 *
 * Ben Prather
 */

#include "decs.hpp"
#include "diffuse.hpp"
#include "self_init.hpp"

#if USE_MPI
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
namespace mpi = boost::mpi;
#endif

#include <cmath>
#include <iostream>
#include <random>
#include <sstream>

using namespace Kokkos;

int main (int argc, char **argv)
{
  size_t ng = 3; // TODO add to dumps/take

#if USE_MPI
  mpi::environment env(argc, argv);
  mpi::communicator world;
#endif
  Kokkos::initialize (argc, argv);
  {
    std::cout << "K/HARM v.alpha" << std::endl;
    std::cout << "Using Kokkos environment:" << std::endl;
    DefaultExecutionSpace::print_configuration(std::cout);
    std::cout << std::endl;

    // TODO make right for parallel HDF5 later
    // HighFive::File input(argv[1], HighFive::File::ReadOnly);
    // auto prims_shape = input.getDataSet("/prims").getDimensions();
    // int n1 = prims_shape[0];
    // int n2 = prims_shape[1];
    // int n3 = prims_shape[2];
    // int nprim = prims_shape[3];
    // std::cout << "Input size: " << n1 << "x" << n2 << "x" << n3 << "x" << nprim << std::endl;

    // Contiguous array of just n1xn2xn3 for input
    // GridPrimsHost h_prims_input("prims_input", n1, n2, n3, nprim);
    // input.getDataSet("/prims").read(h_prims_input.data());

    int n1 = 128, n2 = 128, n3 = 128, nprim = 8;
    GridPrimsHost h_prims_input = mhdmodes(n1, n2, n3, nprim);

    GridPrims prims("prims", n1+2*ng, n2+2*ng, n3+2*ng, nprim);
    GridPrims prims_temp("prims_temp", n1+2*ng, n2+2*ng, n3+2*ng, nprim);
    auto h_prims = Kokkos::create_mirror_view(prims);
    auto h_prims_temp = Kokkos::create_mirror_view(prims);

    MDRangePolicy<OpenMP, Rank<3>> all_range({0,0,0}, {n1,n2,n3});
    Kokkos::parallel_for("diff_all", all_range,
                 KOKKOS_LAMBDA (int i, int j, int k) {
      for (int p=0; p < nprim; ++p) h_prims(i+ng,j+ng,k+ng,p) = h_prims_input(i,j,k,p);
    });

    // copy TO DEVICE
    Kokkos::deep_copy(prims, h_prims);

    for (int iter=0; iter<1000; iter++) {
      diffuse_all(prims, prims_temp);
    }



  }
  Kokkos::finalize(); // TODO additionally call on exceptions?

  return 0;
}
