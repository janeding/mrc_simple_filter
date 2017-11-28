// (Note: For gcc version 4.8.3, you must compile using: g++ -std=c++11)

#include <cassert>
#include <cmath>
#include <iostream>
#include <string>
#include <fftw3.h>
using namespace std;
#include <alloc3d.h>
#include <filter3d.h>
#include <threshold.h>
#include <mrc_simple.h>
#include "settings.h"





// Generate a filter
//
//   h(x,y,z) = h_xy(r) * h_z(z)
//
// In the XY plane, the filter used is:
//
//   h_xy(r) = A*exp(-(r/s)^m) - B*exp(-(r/t)^n)
//
// (Here, "r" refers to disance from the Z axis)
// (and "A" and "B" are chosen so that h_xy(r=0) = A+B = 1, and
//  the sum of h_xy(r) for all pixels in the filter range is zero.)
// Along the Z direction, the filter used is:
//
//   h_z(r) = C*exp(-0.5*(z/u)^2)
//
// (Where "C" is chosen so the (discrete) sum-over-z of  h_z(z) is 1)





template<class RealNum, class Integer>
// Generate a generalized "difference of Gaussians" filter of the form:
//    h_xy(r) = A*exp(-(r/s)^m) - B*exp(-(r/t)^n)
float TryGenFilter2D(Integer size[2],
                     RealNum **aaf,
                     RealNum A,           //"A" coeff in formula
                     RealNum B,           //"B" coeff in formula
                     RealNum width_s,     //"s" parameter in formula
                     RealNum width_t,     //"t" parameter in formula
                     RealNum exponent_m,  //"m" parameter in formula
                     RealNum exponent_n,  //"n" parameter in formula
                     RealNum cutoff = 0)
{
  RealNum total = 0;
  for (Integer iy=-size[1]; iy<=size[1]; iy++) {
    for (Integer ix=-size[0]; ix<=size[0]; ix++) {
      RealNum r = sqrt(ix*ix + iy*iy);
      RealNum h = (A * exp(-pow(r/width_s, exponent_m)) -
                 B * exp(-pow(r/width_t, exponent_n)));
      if (B * exp(-pow(r/width_t, exponent_n)) < cutoff)
        h = 0.0;
      aaf[iy+size[1]][ix+size[0]] = h;
      total += h;
    }
  }
  return total;
} //TryGenFilter2D()



template<class RealNum, class Integer>
float GenFilter2D(Integer size[2],
                  RealNum **aaf,
                  RealNum width_s,     //"s" parameter in formula
                  RealNum width_t,     //"t" parameter in formula
                  RealNum exponent_m,  //"m" parameter in formula
                  RealNum exponent_n,  //"n" parameter in formula
                  RealNum cutoff = 0.0,
                  RealNum tolerance = 1.0e-6,
                  int   max_iters = 100)
{
  RealNum A = 1.0;
  RealNum Bmin = -1.0;
  RealNum Bmax = 0.0;
  RealNum B = Bmax;
  RealNum delta = 1.0;
  int niters = 0;

  // Search for a value of B which causes the area under
  // the curve to be close to 0
  while ((delta > tolerance) &&
         (niters < max_iters)) {
    delta = TryGenFilter2D(size,
                           aaf,
                           A,
                           B,
                           width_s,
                           width_t,
                           exponent_m,
                           exponent_n,
                           cutoff);
    if (delta > 0) {
      Bmin = B;
      B = 0.5 * (Bmin + Bmax);
    }
    else if (delta < 0) {
      Bmax = B;
      B = 0.5 * (Bmin + Bmax);
    }
    niters++;
  } // while (delta > tolerance) ...
  if (niters > max_iters)
    throw string("Error: The size of the mask does not match the size of the tomogram.\n");
  // The filter at the central peak should have a height of 1
  // Instead it has a neight of A + B.
  // Rescale all entries by dividing by A + B  (Keep in mind "B" is negative)
  assert(A > -B);
  for (Integer iy=-size[1]; iy<=size[1]; iy++)
    for (Integer ix=-size[0]; ix<=size[0]; ix++)
      aaf[iy+size[1]][ix+size[0]] /= (A + B);
  A /= (A+B);
  B /= (A+B);
} //GenFilter2D()



template<class RealNum, class Integer>
void GenFilterGauss1D(Integer size,
                      RealNum *af,
                      RealNum sigma,
                      RealNum cutoff = 0.0)
{
  RealNum sum = 0.0;
  for (Integer i=-size; i<=size; i++) {
    RealNum h = exp(-0.5*(i*i)/(sigma*sigma));
    sum += h;
    af[i+size] = h;
  }

  //Normalize:
  for (Integer i=-size; i<=size; i++)
    af[i+size] /= sum;
} //GenFilterGauss1D()





template<class RealNum, class Integer>
void GenFilter3D(Integer size[3],
                 RealNum ***aaaf,
                 RealNum width_s,     //"s" parameter in formula
                 RealNum width_t,     //"t" parameter in formula
                 RealNum width_u,     //"t" parameter in formula
                 RealNum exponent_m,  //"m" parameter in formula
                 RealNum exponent_n,  //"n" parameter in formula
                 RealNum window_sigma_cutoff_xy,
                 RealNum window_sigma_cutoff_z)
{

  // Allocate some temporary tables we will need:
  RealNum **aaf_xy = new RealNum* [size[1] * 2 + 1];
  for (Integer iy=-size[1]; iy<=size[1]; iy++)
    aaf_xy[iy] = new RealNum [size[0]];
  RealNum *af_z = new RealNum [size[2] * 2 + 1];

  // Generate the filter in the xy direction:

  RealNum size_xy[2] = {size[0], size[1]}; //(unnecessary, but makes code more clear)

  // Generate the filter in the Z direction
  GenFilterGauss1D(size[2],
                   af_z,
                   width_u,     //"s" parameter in formula
                   window_sigma_cutoff_z);


  // The 3-D filter is the product of the two filters:

  RealNum total = 0;
  for (Integer iz=-size[2]; iz<=size[2]; iz++) {
    for (Integer iy=-size[1]; iy<=size[1]; iy++) {
      for (Integer ix=-size[0]; ix<=size[0]; ix++) {
        RealNum h = aaf_xy[ix][iy] * af_z[iz];
        total += h;
      }
    }
  }

  // Deallocate temporary tables:
  for (Integer iy=-size[1]; iy<=size[1]; iy++)
    delete [] aaf_xy[iy];
  delete [] aaf_xy;
  delete [] af_z;
} //GenFilter3D()







int main(int argc, char **argv) {
  try {
    Settings settings; // parse the command-line argument list from the shell
    settings.ParseArgs(argc, argv);

    // Read the input tomogram
    cerr << "Reading tomogram \""<<settings.in_file_name<<"\"" << endl;
    MrcSimple tomo;
    tomo.Read(settings.in_file_name, false);
    // (Note: You can also use "tomo.Read(cin);" or "cin >> tomo;")
    tomo.PrintStats(cerr);      //Optional (display the tomogram size & format)

    // ---- mask ----

    // Optional: if there is a "mask", read that too
    MrcSimple mask;
    if (settings.mask_file_name != "") {
      cerr << "Reading mask \""<<settings.mask_file_name<<"\"" << endl;
      mask.Read(settings.mask_file_name, false);
      if ((mask.mrc_header.nvoxels[0] != tomo.mrc_header.nvoxels[0]) ||
          (mask.mrc_header.nvoxels[1] != tomo.mrc_header.nvoxels[1]) ||
          (mask.mrc_header.nvoxels[2] != tomo.mrc_header.nvoxels[2]))
        throw string("Error: The size of the mask does not match the size of the tomogram.\n");
      // The mask should be 1 everywhere we want to consider, and 0 elsewhere.
      if (settings.use_mask_select) {
        for (int iz=0; iz<mask.mrc_header.nvoxels[2]; iz++)
          for (int iy=0; iy<mask.mrc_header.nvoxels[1]; iy++)
            for (int ix=0; ix<mask.mrc_header.nvoxels[0]; ix++)
              if (mask.aaafDensity[iz][iy][ix] == settings.mask_select)
                mask.aaafDensity[iz][iy][ix] = 1.0;
              else
                mask.aaafDensity[iz][iy][ix] = 0.0;
      }
    }

    if (settings.in_rescale01)
      tomo.Rescale01(mask.aaafDensity);

    // ---- make an array that will store the new tomogram we will create ----

    cerr << "allocating space for new tomogram..." << endl;
    MrcSimple out_tomo = tomo; //this will take care of allocating the array

    float voxel_width[3] = {1.0, 1.0, 1.0};

    // ---- filtering ----
    if (settings.voxel_width > 0.0) {
      voxel_width[0] = settings.voxel_width;
      voxel_width[1] = settings.voxel_width;
      voxel_width[2] = settings.voxel_width;
    }
    else {
      voxel_width[0] = 0.1*tomo.mrc_header.cellA[0]/tomo.mrc_header.nvoxels[0];
      voxel_width[1] = 0.1*tomo.mrc_header.cellA[1]/tomo.mrc_header.nvoxels[1];
      voxel_width[2] = 0.1*tomo.mrc_header.cellA[2]/tomo.mrc_header.nvoxels[2];
      cerr << "voxel width in nm = ("
           << voxel_width[0] << ", "
           << voxel_width[1] << ", "
           << voxel_width[2] << ")\n";
    }

    if ((voxel_width[0] <= 0.0) ||
        (voxel_width[1] <= 0.0) ||
        (voxel_width[2] <= 0.0))
      throw string("Error in tomogram header: Invalid voxel width(s).\n"
                   "Use the -w argument to specify the voxel width in nm.");

    if (abs((voxel_width[0] - voxel_width[1])
            /
            (0.5*(voxel_width[0] + voxel_width[1]))) > 0.0001)
      throw string("Error in tomogram header: Unequal voxel widths in the x and y directions.\n"
                   "Use the -w argument to specify the voxel width in nm.");
    settings.width_s /= 0.5*(voxel_width[0]+voxel_width[1]);
    settings.width_t /= 0.5*(voxel_width[0]+voxel_width[1]);
    settings.width_u /= voxel_width[2];
    settings.window_sigma_cutoff_xy /= 0.5*(voxel_width[0]+voxel_width[1]);
    settings.window_sigma_cutoff_z /= voxel_width[2];

    if (settings.window_sigma_cutoff_xy > 0.0) {
      settings.window_size[0] = ceil(settings.window_sigma_cutoff_xy * 
				     settings.width_t);
      settings.window_size[1] = ceil(settings.window_sigma_cutoff_xy * 
				     settings.width_t);
    }

    if (settings.window_sigma_cutoff_z > 0.0) {
      settings.window_size[2] = ceil(settings.window_sigma_cutoff_z * 
				     settings.width_u);
    }

    cerr << "applying filter (window size "
	 << settings.window_size[0] << ","
	 << settings.window_size[1] << ","
	 << settings.window_size[2] << ")"
	 << " ..." << endl;

    Filter3D<float, int> filter(settings.window_size);

    GenFilter3D(filter.size,
                filter.aaafWeights,
                settings.width_s,     //"s" parameter in formula
                settings.width_t,     //"t" parameter in formula
                settings.width_u,     //"t" parameter in formula
                settings.exponent_m,  //"m" parameter in formula
                settings.exponent_n,  //"n" parameter in formula
                settings.window_sigma_cutoff_xy,
                settings.window_sigma_cutoff_z);

    filter.Apply(tomo.mrc_header.nvoxels,
                 tomo.aaafDensity,
                 out_tomo.aaafDensity,
                 mask.aaafDensity,
                 true,
                 &cout);

    if (settings.in_rescale01)
      out_tomo.Rescale01(mask.aaafDensity);

    // ----- thresholding and masking: -----
    
    if (settings.use_thresholds) {

      cerr << "thresholding..." << endl;

      for (int iz=0; iz<out_tomo.mrc_header.nvoxels[2]; iz++) {
        for (int iy=0; iy<out_tomo.mrc_header.nvoxels[1]; iy++) {
          for (int ix=0; ix<out_tomo.mrc_header.nvoxels[0]; ix++) {
            if (! settings.use_dual_thresholds)
              out_tomo.aaafDensity[iz][iy][ix] =
                Threshold2(out_tomo.aaafDensity[iz][iy][ix],
                           settings.out_threshold_01_a,
                           settings.out_threshold_01_b);
            else
              out_tomo.aaafDensity[iz][iy][ix] =
                Threshold4(out_tomo.aaafDensity[iz][iy][ix],
                           settings.out_threshold_01_a,
                           settings.out_threshold_01_b,
                           settings.out_threshold_10_a,
                           settings.out_threshold_10_b);
          }
        }
      }
    }

    // Also, after thresholding, check again to see if the mask is zero
    // at this location.  If so, make sure that aaafDensity to 0 there.
    if ((mask.aaafDensity) && (settings.use_mask_out))
      for (int iz=0; iz<mask.mrc_header.nvoxels[2]; iz++)
        for (int iy=0; iy<mask.mrc_header.nvoxels[1]; iy++)
          for (int ix=0; ix<mask.mrc_header.nvoxels[0]; ix++)
            if (mask.aaafDensity[iz][iy][ix] == 0.0)
              out_tomo.aaafDensity[iz][iy][ix] = settings.mask_out;


    // Write the file
    if (settings.out_file_name != "") {
      cerr << "writing tomogram (in float mode)" << endl;
      out_tomo.Write(settings.out_file_name);
      //(You can also use "out_tomo.Write(cout);" or "cout<<out_tomo;")
    }
  }

  catch (string s) {
    cerr << s << endl; // In case of file format error, print message and exit
    exit(1);
  }
}







