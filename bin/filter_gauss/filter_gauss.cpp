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


// (Note: For gcc version 4.8.3, you must compile using: g++ -std=c++11)

template<class RealNum, class Integer>
// Generate a filter with Gaussian weights.
// Must specify the width of the Gaussian in all 3 directions (sigma).
void FillGaussian(Integer size[3],
                  RealNum ***aaaf,
                  RealNum sigma[3], 
		  RealNum cutoff = 0)
{
  RealNum total = 0;
  for (Integer iz=-size[2]; iz<=size[2]; iz++) {
    for (Integer iy=-size[1]; iy<=size[1]; iy++) {
      for (Integer ix=-size[0]; ix<=size[0]; ix++) {
        RealNum h = exp(-0.5*((ix*ix)/(sigma[0]*sigma[0]) + 
                            (iy*iy)/(sigma[1]*sigma[1]) + 
                            (iz*iz)/(sigma[2]*sigma[2])));
	if (h < cutoff)
	  h = 0.0;
        aaaf[iz+size[2]][iy+size[1]][ix+size[0]] = h;
        total += h;
      }
    }
  }
  // Normalize: Make sure the sum of the weights is 1.
  for (Integer iz=-size[2]; iz<=size[2]; iz++)
    for (Integer iy=-size[1]; iy<=size[1]; iy++)
      for (Integer ix=-size[0]; ix<=size[0]; ix++)
        aaaf[iz+size[2]][iy+size[1]][ix+size[0]] /= total;
}



template<class RealNum, class Integer>
// Probably unnecessary:
// Generate a filter with Gaussian weights using a slighly more general formula.
// (In this version, the Gaussian need not be oriented along the x,y,z axis.
//  Instead, the "sigma" parameter can be a matrix with off-diagonal elements.)
void FillGaussian(Integer size[3],
                  RealNum ***aaaf,
                  RealNum sigma[3][3],
		  RealNum cutoff = 0)
{
  RealNum total = 0;
  Integer iLoc[3];
  for (iLoc[2]=-size[2]; iLoc[2]<=size[2]; iLoc[2]++) {
    for (iLoc[1]=-size[1]; iLoc[1]<=size[1]; iLoc[1]++) {
      for (iLoc[0]=-size[0]; iLoc[0]<=size[0]; iLoc[0]++) {
        RealNum sum = 0.0;
        for (Integer i = 0; i < 3; i++)
          for (Integer j = 0; j < 3; j++)
            sum += (iLoc[i]*iLoc[j])/(sigma[i]*sigma[j]);
        RealNum h = exp(-0.5*sum);
	if (h < cutoff)
	  h = 0.0;
        aaaf[iLoc[2]+size[2]]
            [iLoc[1]+size[1]]
            [iLoc[0]+size[0]] = h;
        total += h;
      }
    }
  }
  // Normalize: Make sure the sum of the weights is 1.
  for (Integer iz=-size[2]; iz<=size[2]; iz++)
    for (Integer iy=-size[1]; iy<=size[1]; iy++)
      for (Integer ix=-size[0]; ix<=size[0]; ix++)
        aaaf[iz+size[2]][iy+size[1]][ix+size[0]] /= total;
}






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

    // ---- filtering ----

    if (settings.window_sigma_nm[0] > 0) {
      settings.window_sigma[0] = settings.window_sigma_nm[0] / (0.1*tomo.mrc_header.cellA[0]/tomo.mrc_header.nvoxels[0]);
      settings.window_sigma[1] = settings.window_sigma_nm[1] / (0.1*tomo.mrc_header.cellA[1]/tomo.mrc_header.nvoxels[1]);
      settings.window_sigma[2] = settings.window_sigma_nm[2] / (0.1*tomo.mrc_header.cellA[2]/tomo.mrc_header.nvoxels[2]);
    }
    if (settings.window_sigma_cutoff > 0.0) {
      settings.window_size[0] = ceil(settings.window_sigma_cutoff * 
				     settings.window_sigma[0]);
      settings.window_size[1] = ceil(settings.window_sigma_cutoff * 
				     settings.window_sigma[1]);
      settings.window_size[2] = ceil(settings.window_sigma_cutoff * 
				     settings.window_sigma[2]);
    }

    cerr << "applying filter (window size "
	 << settings.window_size[0] << ","
	 << settings.window_size[1] << ","
	 << settings.window_size[2] << ")"
	 << " ..." << endl;

    Filter3D<float, int> filter(settings.window_size);
    // In this version of the program, assume the user wants a Gaussian filter
    FillGaussian(filter.size,
		 filter.aaafWeights,
		 settings.window_sigma,
		 settings.window_sigma_cutoff_exp);

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







