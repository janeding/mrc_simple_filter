// (Note: For gcc version 4.8.3, you must compile using: g++ -std=c++11)

#include <cassert>
#include <cmath>
#include <iostream>
#include <string>
//#include <fftw3.h>  not needed yet
using namespace std;
#include <alloc.h>
#include <filter.h>
#include <threshold.h>
#include <mrc_simple.h>
#include "settings.h"





template<class RealNum>
inline RealNum SQR(RealNum x) { return x*x; }

template<class RealNum>
inline RealNum MAX(RealNum x, RealNum y) { return ((x<y) ? y : x); }

template<class RealNum, class Integer>




template<class RealNum, class Integer>
// Fill the 2-D array "aaf" with a "generalized Gaussian" function of the form:
//    h_xy(r) = A*exp(-r^m)
// where   r  = sqrt((x/s_x)^2 + (y/s_y)^2)
//   and   A  is determined by normalization of the discrete sum
// (The function returns "A" to the caller.)
RealNum GenFilterGenGauss2D(Integer const size[2],
                            RealNum **aaf,
                            RealNum width[2],    //"s_x", "s_y" parameters
                            RealNum exponent_m,  //"m" parameter in formula
                            RealNum exponent_n,  //"n" parameter in formula
                            RealNum cutoff = 0)
{
  RealNum total = 0;
  for (Integer iy=-size[1]; iy<=size[1]; iy++) {
    for (Integer ix=-size[0]; ix<=size[0]; ix++) {
      RealNum r = sqrt(SQR(ix/width[0]) + SQR(iy/width[1]));
      RealNum h = exp(-pow(r, exponent_m));
      if (h < cutoff)
        h = 0.0;
      aaf[iy+size[1]][ix+size[0]] = h;
      total += h;
    }
  }
  for (Integer iy=-size[1]; iy<=size[1]; iy++) {
    for (Integer ix=-size[0]; ix<=size[0]; ix++) {
      aaf[iy+size[1]][ix+size[0]] /= total;
    }
  }
  return 1.0 / total;
} //GenFilterGenGauss2D()



template<class RealNum, class Integer>

// Fill the "aaf" array with a difference of (generalized) Gaussians:
RealNum GenFilterGenDog2D(Integer const size[2], //Size of the filter we want
                          RealNum **aaf,       //array of filter values (h(x,y))
                          RealNum width_a[2],  //"a" parameter in formula
                          RealNum width_b[2],  //"b" parameter in formula
                          RealNum exponent_m,  //"m" parameter in formula
                          RealNum exponent_n,  //"n" parameter in formula
                          RealNum tolerance = 1.0e-6,
                          int     max_iters = 100)
{
  RealNum A = 1.0;
  RealNum Bmin = -1.0;
  RealNum Bmax = 0.0;
  RealNum B = Bmax;
  RealNum delta = 1.0;
  int niters = 0;
  RealNum **aafA, **aafB;

  // Optional:
  // Set the filter to zero whenever the value decays below "cutoff"
  // Make sure "cutoff" is compatible with the cutoffs in the x,y directions
  // This gives the filter a round shape (instead of a rectangular shape).
  RealNum cutoff = 1.0;
  for (int d=0; d<2; d++) {
    RealNum h;
    h = exp(-pow(window_cutoff_ratio[d], exponent_m));
    if (h < cutoff)
      cutoff = h;

    h = exp(-pow(window_cutoff_ratio[d], exponent_n));
    if (h < cutoff)
      cutoff = h;
  }

  filterXY_A = Filter2D<float, int>(filter_window_size);
  filterXY_B = Filter2D<float, int>(filter_window_size);
  A = GenFilterGenGauss2D(settings.window_size,
                          filterXY_A.aaf,
                          settings.width_a,    //"a_x", "a_y" parameters
                          settings.exponent_m, //"m" parameter in formula
                          cutoff);
  B = GenFilterGenGauss2D(settings.window_size,
                          filterXY_B.aaf,
                          settings.width_b,    //"b_x", "b_y" parameters
                          settings.exponent_n, //"n" parameter in formula
                          cutoff);

  // The "difference of gaussians" filter is the difference between
  // these two (generalized) gaussian filters.

  // (We also rescale the filter upwards by 1.0/(A-B) to insure that the central
  //  peak has height 1.  (The total area under the filter will be 0)).
  for (Integer iy=-size[1]; iy<=size[1]; iy++) {
    for (Integer ix=-size[0]; ix<=size[0]; ix++) {
      aaf[iy+size[1]][ix+size[0]] = (
                                     (filterXY_A.aaf[iy+size[1]][ix+size[0]]
                                      -
                                      filterXY_B.aaf[iy+size[1]][ix+size[0]])
                                     /
                                     (A - B)
                                     );
    }
  }
} //GenFilterGenDog2D()




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

void
ApplyGauss3D(Integer const filter_window_size[3],
             Integer const image_size[3], 
             RealNum ***aaafSource,
             RealNum ***aaafDest,
             RealNum ***aaafMask,
             bool precompute_mask_times_source = false)
{
  //Allocate filters in all 3 directions.  (Later apply them sequentially.)
  aFilter = new Filter1D<float, int> [3];
  for (int d=0; d <= 3; d++) {
    aFilter[d].Resize(filter_window_size[d]);
    GenFilterGauss1D(aFilter[d].size,
                     aFilter[d].afWeights,
                     widths[d]);
  }

  // Create temporary 1-D arrays to perform the filter in each direction:
  float **aafSource[3];
  float **aafDest[3];
  float **aafMask[3];
  for (int d=0; d <= 3; d++) {
    aafDest[d] = new float image_size[d];
    aafSource[d] = new float image_size[d];
    aafMask[d]  = NULL;
    if (aaafMask)
      aafMask[d] = new float image_size[d];
  }

  // Initialize the aaafDest array
  for (Integer iz = 0; iz < image_size[2]; iz++)
    for (Integer iy = 0; iy < image_size[1]; iy++)
      for (Integer ix = 0; ix < image_size[0]; ix++)
        aaafDest[iz][iy][ix] = aaafSource[iz][iy][ix];

  // This is a "separable" filter.
  // You can apply the filter sequentially, in the X, Y, Z directions
  // (instead of applying the filter simultaneously in all 3 directions,
  //  which would be much slower)

  // Apply the filter in the Z direction (d=2):
  int d=2;
  for (Integer iy = 0; iy < image_size[1]; iy++) {
    for (Integer ix = 0; ix < image_size[0]; ix++) {
      for (Integer iz = 0; iz < image_size[2]; iz++) {
        aafSource[d][iz] = aaafDest[iz][iy][ix];
        if (afMask)
          aafMask[d][iz] = aaafMask[iz][iy][ix];
      }
      aFilter[d].Apply(image_size[d],
                       aafSource[d],
                       aafDest[d],
                       aafMask[d],
                       precompute_mask_times_source);
      for (Integer iz = 0; iz < image_size[d]; iz++) {
        //It would be wasteful to allocate a temporary array to store this
        //Instead store the result of the convolution in the original array:
        aaafDest[iz][iy][ix] = aafDest[d][iz];
      }
    } //for (Integer ix = 0; ix < image_size[0]; ix++)
  } //for (Integer iy = 0; iy < image_size[1]; iy++)


  // Apply the filter in the Y direction:
  int d=1;
  for (Integer iz = 0; iz < image_size[2]; iz++) {
    for (Integer ix = 0; ix < image_size[0]; ix++) {
      for (Integer iy = 0; iy < image_size[1]; iy++) {
        aafSource[d][iy] = aaafDest[iz][iy][ix];
        if (afMask)
          aafMask[d][iy] = aaafMask[iz][iy][ix];
      }
      aFilter[d].Apply(image_size[d],
                       aafSource[d],
                       aafDest[d],
                       aafMask[d],
                       precompute_mask_times_source);
      for (Integer iy = 0; iy < image_size[d]; iy++) {
        //It would be wasteful to allocate a temporary array to store this
        //Instead store the result of the convolution in the original array:
        aaafDest[iz][iy][ix] = aafDest[d][iy];
      }
    } //for (Integer ix = 0; ix < image_size[0]; ix++) {
  } //for (Integer iz = 0; iz < image_size[2]; iz++)

  // Apply the filter in the X direction:
  int d=0;
  for (Integer iz = 0; iz < image_size[2]; iz++) {
    for (Integer iy = 0; iy < image_size[1]; iy++) {
      for (Integer ix = 0; ix < image_size[0]; ix++) {
        aafSource[d][ix] = aaafDest[iz][iy][ix];
        if (afMask)
          aafMask[d][ix] = aaafMask[iz][iy][ix];
      }
      aFilter[d].Apply(image_size[d],
                       aafSource[d],
                       aafDest[d],
                       aafMask[d],
                       precompute_mask_times_source);
      for (Integer ix = 0; ix < image_size[d]; ix++) {
        //It would be wasteful to allocate a temporary array to store this
        //Instead store the result of the convolution in the original array:
        aaafDest[iz][iy][ix] = aafDest[d][ix];
      }
    } //for (Integer iy = 0; iy < image_size[1]; iy++)
  } //for (Integer iz = 0; iz < image_size[2]; iz++)


  for (int d=0; d<3; d++) {
    delete [] aafSource[d];
    delete [] aafDest[d];
    if (aafMask[d])
      delete [] aafMask[d];
  }
  delete [] aFilter;
} //ApplyGauss3D()






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
      // Did the user manually specify the width of each voxel?
      voxel_width[0] = settings.voxel_width;
      voxel_width[1] = settings.voxel_width;
      voxel_width[2] = settings.voxel_width;
    }
    else {
      // Otherwise, infer it from the header of the MRC file
      voxel_width[0] = tomo.mrc_header.cellA[0]/tomo.mrc_header.nvoxels[0];
      voxel_width[1] = tomo.mrc_header.cellA[1]/tomo.mrc_header.nvoxels[1];
      voxel_width[2] = tomo.mrc_header.cellA[2]/tomo.mrc_header.nvoxels[2];
      if (settings.voxel_width_divide_by_10) {
        voxel_width[0] *= 0.1;
        voxel_width[1] *= 0.1;
        voxel_width[2] *= 0.1;
      }
      cerr << "voxel width in physical units = ("
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
    for (int d=0; d<3; d++) {
      settings.width_a[d] /= voxel_width[d];
      settings.width_b[d] /= voxel_width[d];
      //settings.window_cutoff_ratio[d] /= voxel_width[d];
      //settings.window_cutoff_dog[d] /= voxel_width[d];
      //settings.window_cutoff_gauss[d] /= voxel_width[d];
      settings.filter_window_size[d] = ceil(settings.window_cutoff_ratio[d] * 
                                            MAX(settings.width_a[d],
                                                settings.width_b[d]));
    }

    cerr << "applying filter (window size in voxels: "
         << settings.filter_window_size[0] << ","
         << settings.filter_window_size[1] << ","
         << settings.filter_window_size[2] << ")"
         << " ..." << endl;

 
    if (settings.filter_type = GAUSS) {

      cerr << " Filter Used:\n"
        " h(x,y,z)   = A*exp(-0.5*((x/s_x)^2 + (y/s_y)^2 + (z/s_z)^))\n"
        " ... where  A,  B,  a_x,  b_x  equal:\n"
           << " " << A << " " << B
           << " " << width_a[0]/sqrt(2)
           << " " << width_b[0]/sqrt(2) << "\n";
        " ...   and A,  B,  a_y,  b_y  equal:\n"
           << " " << A << " " << B
           << " " << width_a[1]/sqrt(2)
           << " " << width_b[1]/sqrt(2) << "\n";
        " ...   and A,  B,  a_z,  b_z  equal:\n"
           << " " << A << " " << B
           << " " << width_a[2]/sqrt(2)
           << " " << width_b[2]/sqrt(2) << "\n"
          " You can plot this function in the X,Y, or Z directions using:\n"
          " draw_filter_1D.py -dog  A  B  a  b\n";
        
      ApplyGauss3D(settings.filter_window_size,
                   tomo.mrc_header.nvoxels,
                   tomo.aaafDensity,
                   out_tomo.aaafDensity,
                   mask.aaafDensity,
                   true);

    } //if (settings.filter_type = GAUSS)

    if (settings.filter_type = DOG_FAST) {

      float ***aaafTemp; //temporary array to store partially processed tomogram
      float *afTemp;     //temporary array to store partially processed tomogram
      int   *size = tomo.mrc_header.nvoxels;
      Alloc3D(tomo.mrc_header.nvoxels,
              &afTemp,
              &aaafTemp);

      // Initialize the aaafTemp array
      for (int iz = 0; iz < size[2]; iz++)
        for (int iy = 0; iy < size[1]; iy++)
          for (int ix = 0; ix < size[0]; ix++)
            aaafTemp[iz][iy][ix] = aaafSource[iz][iy][ix];

      cerr << " Filter Used:\n"
        " h(x,y,z)   = h_a(x,y,z) - h_b(x,y,z)\n"
        " h_a(x,y,z) = A*exp(-0.5*((x/a_x)^2 + (y/a_y)^2 + (z/a_z)^2))\n"
        " h_b(x,y,z) = B*exp(-0.5*((x/b_x)^2 + (y/b_y)^2 + (z/b_z)^2))\n"
        " ... where  A,  B,  a_x,  b_x  equal:\n"
           << " " << A << " " << B
           << " " << width_a[0]/sqrt(2)
           << " " << width_b[0]/sqrt(2) << "\n";
        " ...   and A,  B,  a_y,  b_y  equal:\n"
           << " " << A << " " << B
           << " " << width_a[1]/sqrt(2)
           << " " << width_b[1]/sqrt(2) << "\n";
        " ...   and A,  B,  a_z,  b_z  equal:\n"
           << " " << A << " " << B
           << " " << width_a[2]/sqrt(2)
           << " " << width_b[2]/sqrt(2) << "\n"
          " You can plot this function in the X,Y, or Z directions using:\n"
          " draw_dog_filter_1D.py  A  B  a  b\n";
        
      // Rather than 
      // Convolve the original source with the 1st Gaussian
      ApplyGauss3D(tomo.mrc_header.nvoxels,
                   tomo.aaafDensity,
                   out_tomo.aaafDensity,
                   mask.aaafDensity,
                   true);

      // Convolve the original source with the 2nd Gaussian
      ApplyGauss3D(tomo.mrc_header.nvoxels,
                   tomo.aaafDensity,
                   aaafTemp,
                   mask.aaafDensity,
                   true);

      // Subtract the second convolved signal from the first
      for (int iz = 0; iz < size[2]; iz++)
        for (int iy = 0; iy < size[1]; iy++)
          for (int ix = 0; ix < size[0]; ix++)
            out_tomo.aaafDensity[iz][iy][ix] -= aaafTemp[iz][iy][ix];

      // Deallocate the temporary array
      Dealloc3D(tomo.mrc_header.nvoxels,
                &afTemp,
                &aaafTemp);

    } //if (settings.filter_type = DOG_FAST)

    else if (settings.filter_type = DOG) {
      // Optional:
      // Set the filter to zero whenever the value decays below "cutoff"
      // Make sure "cutoff" is compatible with the cutoffs in the x,y directions
      // This gives the filter a round shape (instead of a rectangular shape).
      cutoff = 1.0;
      for (int d=0; d<3; d++) {
        if (settings.width_a[d] > settings.width_b[d])
          h = exp(-pow(settings.window_cutoff_ratio[d], settings.exponent_m));
        else
          h = exp(-pow(settings.window_cutoff_ratio[d], settings.exponent_n));
        if (h < cutoff)
          cutoff = h;
      }
      throw string("Error: The general 3-D DOG filter (supporting exponents m,n!=2)\n"
                   "       has not been implemented yet.\n"
                   "      (However implementing should be trivial to do.\n"
                   "       Edit the code to define a functions \"GenFilterGenDog3D()\" and\n"
                   "       GenFilterGenGause3D This should be easy. They will be nearly identical\n"
                   "       to the functions \"GenFilterGenDog2D()\" and \"GenFilterGenGauss2D()\")\n"
                   "For now, you can use an ordinary 3-D DOG filter with default exponents m=n=2.\n");

      Filter3D<float, int> filter(settings.filter_window_size);
      //GenFilterGenDog3D(filter.size,
      //                  filter.aaafWeights,
      //                  settings.width_a,     //"a" parameter in formula
      //                  settings.width_b,     //"b" parameter in formula
      //                  settings.exponent_m,  //"m" parameter in formula
      //                  settings.exponent_n,  //"n" parameter in formula
      //                  settings.window_cutoff_ratio_dog);
    } //else if (settings.filter_type = DOG)


    else if (settings.filter_type = DOGXY) {
      // Generate a filter
      //
      //   h(x,y,z) = h_xy(x,y) * h_z(z)
      //
      // In the XY plane, the filter used is:
      //
      //   h_xy(x,y) = A*exp(-ra^m) - B*exp(-rb^n), where:
      //          ra = sqrt((x/a_x)^2 + (y/a_y)^2)
      //          rb = sqrt((x/b_x)^2 + (y/b_y)^2)
      //
      // (Here, "r" refers to disance from the Z axis)
      // (and "A" and "B" are chosen so that h_xy(r=0) = A+B = 1, and
      //  the sum of h_xy(r) for all pixels in the filter range is zero.)
      // Along the Z direction, the filter used is:
      //
      //   h_z(r) = C*exp(-(1/2)*(z/a_z)^2)
      //
      // (Where "C" is chosen so the (discrete) sum-over-z of  h_z(z) is 1)
      //
      //Take advantage of the fact that the filter we
      //are using (ie. the function we are convolving with the source image),
      //is the product of a function of X,Y, with a function of Z.
      //This makes it a "seperable" filter:  We can perform the filters in the Z
      //direction, followed by filtering the result in the X & Y directions.
      //This reduces the computation by a factor of O(filter.size[2]))
      //(A 1-D convolution followed by a 2-D convolution is much faster per 
      // voxel than a full 3-D convolution.)

      // First, generate the filter in the Z direction:

      filterZ = Filter1D<float, int>(settings.filter_window_size[2]);
      GenFilterGauss1D(filterZ.size,
                       filterZ.afWeights,
                       width_a[2]);



      // Then generate the filter in the XY directions

      float filter_window_size_xy[2] = {settings.filter_window_size[0],  //(unnecessary, but 
                                 settings.filter_window_size[1]}; // makes more clear)
      filterXY = Filter2D<float, int>(filter_window_size_xy);

      GenFilterDog2D(filterXY.size,
                     filterXY.aafWeights,
                     width_a,     //"a" parameter in formula
                     width_b,     //"b" parameter in formula
                     exponent_m,  //"m" parameter in formula
                     exponent_n,  //"n" parameter in formula
                     1.0e-6,
                     100);
  

      // Precompute the effect of the filter in the Z direction.

      // Create temporary 1-D arrays to perform the filter in the Z-direction:
      afDensityZorig = new float tomo.mrc_header.nvoxels[2];
      afDensityZnew  = new float tomo.mrc_header.nvoxels[2];
      afMask         = NULL;
      if (mask.aaafDensity)
        afMask       = new float tomo.mrc_header.nvoxels[2];

      // Then apply the filter in the Z direction
      // (and store the filtered 3-D image in the original tomogram array)
      for (int ix = 0; ix < tomo.mrc_header.nvoxels[0]; ix++) {
        for (int iy = 0; iy < tomo.mrc_header.nvoxels[1]; iy++) {
          for (int iz = 0; iz < tomo.mrc_header.nvoxels[2]; iz++) {
            afDensityZorig[iz] = tomo.aaafDensity[iz][iy][ix];
            if (afMask)
              afMask[iz] = mask.aaafDensity[iz][iy][ix];
          }
          afDensityZorig[iz] = tomo.aaafDensity[iz][iy][ix];
          filterZ.Apply(tomo.mrc_header.nvoxels[2],
                        afDensityZorig,
                        afDensityZnew,
                        afMask,
                        true);
          //It would be wasteful to allocate a temporary array to store this
          //Instead store the result of the convolution in the original array:
          for (int iz = 0; iz < tomo.mrc_header.nvoxels[2]; iz++)
            tomo.aaafDensity[iz][iy][ix] = afDensityZnew[iz];
        } //for (int iy = 0; iy < tomo.mrc_header.nvoxels[1]; iy++) {
      } // for (int ix = 0; ix < tomo.mrc_header.nvoxels[0]; ix++) {
      delete [] afDensityZorig;
      delete [] afDensityZnew;
      if (afMask)
        delete [] afMask;

      // Now apply the filter in the X and Y directions:
      for (int iz = 0; iz < tomo.mrc_header.nvoxels[2]; iz++) {
        filterXY.Apply(tomo.mrc_header.nvoxels,
                       tomo.aaafDensity[iz],
                       out_tomo.aaafDensity[iz],
                       mask.aaafDensity[iz],
                       true,
                       &cout);
      }
    } //else if (settings.filter_type = DOGXY)






    // --- Rescale so that the lowest, highest voxels have density 0 and 1? ---

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







