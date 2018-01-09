// (Note: For gcc version 4.8.3, you must compile using: g++ -std=c++11)

#include <cassert>
#include <cmath>
#include <iostream>
#include <string>
//#include <fftw3.h>  not needed yet
using namespace std;
#include <err_report.h>
#include <alloc2d.h>
#include <alloc3d.h>
#include <filter1d.h>
#include <filter2d.h>
#include <filter3d.h>
#include <threshold.h>
#include <mrc_simple.h>
#include <err_report.h>
#include "settings.h"




template<class RealNum >
inline RealNum ABS(RealNum x) { return ((x<0.0) ? -x: x); }

template<class RealNum >
inline RealNum SQR(RealNum x) { return x*x; }

template<class RealNum >
inline RealNum MAX(RealNum x, RealNum y) { return ((x<y) ? y : x); }


template<class RealNum>
// GenFilterGauss1D generates a 1-D filter and fills its array with values
// corresponding to a normalized Gaussian evaluated at evenly spaced intervals.
// The caller must specify the "sigma" parameter (width of the Gaussian,
// in units of pixels/voxels), in addition to the "halfwidth" parameter, which
// indicates the number of entries in the array (in units of pixels/voxels).
Filter1D<RealNum, int>
GenFilterGauss1D(RealNum sigma,     // The "sigma" paramgeter in the Gaussian
                 int halfwidth) // number of entries in the filter array / 2
{
  Filter1D<RealNum, int> filter(halfwidth);

  RealNum sum = 0.0;
  for (int i=-halfwidth; i<=halfwidth; i++) {
    if (sigma == 0.0) //(don't crash when sigma=0)
      filter.afWeights[i+halfwidth] = ((i == 0) ? 1.0 : 0.0);
    else
      filter.afWeights[i+halfwidth] = exp(-0.5*(i*i)/(sigma*sigma));
    sum += filter.afWeights[i+halfwidth];
  }

  //Normalize:
  for (int i=-halfwidth; i<=halfwidth; i++) {
    filter.afWeights[i+halfwidth] /= sum;
    cerr <<"afWeights["<<i<<"] = "  //FOR DEBUGGING REMOVE EVENTUALLY
         << filter.afWeights[i+halfwidth] << endl;
  }
  return filter;
} //GenFilterGauss1D(sigma, halfwidth)


template<class RealNum>
// This function generates a 1-D filter and fills its array with values
// corresponding to a normalized Gaussian evaluated at evenly spaced intervals.
// The caller must specify the "sigma" parameter (width of the Gaussian).
// and a fractional number between 0 and 1 which indicate how far the
// filter can decay.  Only voxels whose Gaussian intensity decays by less
// than this threshold (relative to the central peak) will be kept.
Filter1D<RealNum, int>
GenFilterGauss1DThresh(RealNum sigma,
                       RealNum window_threshold)
{
  // How wide should the filter be?
  int halfwidth;
  assert(window_threshold > 0.0);
  // Choose the filter domain window based on the "filter_cutoff_threshold"
  //    window_threshold = exp(-0.5*(halfwidth/sigma)^2);
  //    -> (halfwidth/sigma)^2 = -2*log(window_threshold)
  halfwidth = floor(sigma * sqrt(-2*log(window_threshold)));

  return GenFilterGauss1D(sigma, halfwidth);
} //GenFilterGauss1D(sigma, window_threshold)



template<class RealNum>
// Create a 2-D filter and fill it with a "generalized Gaussian" function:
//    h_xy(r) = A*exp(-r^m)
// where   r  = sqrt((x/s_x)^2 + (y/s_y)^2)
//   and   A  is determined by normalization of the discrete sum
// Note: "A" is equal to the value stored in the middle of the array,
//       The caller can determine what "A" is by looking at this value.
Filter2D<RealNum, int>
GenFilterGenGauss2D(RealNum width[2],    //"s_x", "s_y" parameters
                    RealNum m_exp,       //"m" exponent parameter
                    int halfwidth[2])
{
  RealNum window_threshold = 1.0;
  for (int d=0; d<2; d++) {
    RealNum h = ((width[d]>0) ? exp(-pow(halfwidth[d]/width[d], m_exp)) : 1.0);
    if (h < window_threshold)
      window_threshold = h;
  }
  Filter2D<RealNum, int> filter(halfwidth);
  RealNum total = 0;
  for (int iy=-halfwidth[1]; iy<=halfwidth[1]; iy++) {
    for (int ix=-halfwidth[0]; ix<=halfwidth[0]; ix++) {
      RealNum r = sqrt(SQR(ix/width[0]) + SQR(iy/width[1]));
      RealNum h = ((r>0) ? exp(-pow(r, m_exp)) : 1.0);
      if (ABS(h) < window_threshold)
        h = 0.0; // this eliminates corner entries which fall below threshold
                 // (and eliminates anisotropic artifacts due to these corners)
                 // There's no reason to keep any entries less than min value.
      filter.aafWeights[iy+halfwidth[1]][ix+halfwidth[0]] = h;
      total += h;
    }
  }
  // normalize:
  for (int iy=-halfwidth[1]; iy<=halfwidth[1]; iy++) {
    for (int ix=-halfwidth[0]; ix<=halfwidth[0]; ix++) {
      filter.aafWeights[iy+halfwidth[1]][ix+halfwidth[0]] /= total;
      //FOR DEBUGGING REMOVE EVENTUALLY
      cerr << "threshold=" << window_threshold
           <<", aafWeights["<<iy<<"]["<<ix<<"] = "
           << filter.aafWeights[iy+halfwidth[1]][ix+halfwidth[0]] << endl;
    }
  }
  return filter;
} //GenFilterGenGauss2D(width, m_exp, halfwidth, window_thresh)


template<class RealNum>
Filter2D<RealNum, int>
GenFilterGenGauss2DThresh(RealNum width[2],    //"s_x", "s_y" parameters
                          RealNum m_exp,       //"m" parameter in formula
                          RealNum window_threshold) //controls window width
{
  // choose the filter window width based on the window_threshold
  int halfwidth[2];
  int ix = 0;

  for (int d=0; d<2; d++) {
    // Choose the filter domain window based on the "filter_cutoff_threshold"
    //    window_threshold = exp(-(halfwidth/sigma)^m_exp);
    //    -> (halfwidth/sigma)^m_exp = -log(window_threshold)
    halfwidth[d] = floor(width[d] * pow(-log(window_threshold), 1.0/m_exp));
  }
  return GenFilterGenGauss2D(width,
                             m_exp,
                             halfwidth);
} //GenFilterGenGauss2DThresh(width, m_exp, window_thresh)


template<class RealNum>
// Create a 2-D filter and fill it with a difference of (generalized) Gaussians:
// This version requires that the caller has already created individual
// filters for the two gaussians.
// All this function does is subtract one filter from the other (and rescale).
Filter2D<RealNum, int> 
_GenFilterGenDog2D(RealNum width_a[2],  //"a" parameter in formula
                   RealNum width_b[2],  //"b" parameter in formula
                   RealNum m_exp,  //"m" parameter in formula
                   RealNum n_exp,  //"n" parameter in formula
                   Filter2D<RealNum, int>& filterXY_A, //filters for the two
                   Filter2D<RealNum, int>& filterXY_B, //gaussians
                   RealNum *pA=NULL, //optional:report A,B coeffs to user
                   RealNum *pB=NULL) //optional:report A,B coeffs to user
{
  RealNum A, B;
  //A, B = height of the central peak
  //       The central peak is located in the middle of the filter's array
  //       (at position "halfwidth")
  A = filterXY_A.aafWeights[filterXY_A.halfwidth[0]][filterXY_A.halfwidth[1]];
  B = filterXY_B.aafWeights[filterXY_B.halfwidth[0]][filterXY_B.halfwidth[1]];


  // The "difference of gaussians" filter is the difference between
  // these two (generalized) gaussian filters.
  int halfwidth[2];
  halfwidth[0] = MAX(filterXY_A.halfwidth[0], filterXY_B.halfwidth[0]);
  halfwidth[1] = MAX(filterXY_A.halfwidth[1], filterXY_B.halfwidth[1]);
  Filter2D<RealNum, int> filter(halfwidth);

  cerr << "Array of 2D filter entries:" << endl;
  for (int iy=-halfwidth[1]; iy<=halfwidth[1]; iy++) {
    for (int ix=-halfwidth[0]; ix<=halfwidth[0]; ix++) {
      filter.aafWeights[iy+halfwidth[1]][ix+halfwidth[0]] = 0.0;

      // The two filters may have different widths, so we have to check
      // that ix and iy lie within the domain of these two filters before
      // adding or subtracting their values from the final DOG filter.
      if (((-filterXY_A.halfwidth[0]<=ix) && (ix<=filterXY_A.halfwidth[0])) &&
          ((-filterXY_A.halfwidth[1]<=iy) && (iy<=filterXY_A.halfwidth[1])))

        filter.aafWeights[iy+halfwidth[1]][ix+halfwidth[0]] +=
          filterXY_A.aafWeights[iy+filterXY_A.halfwidth[1]][ix+filterXY_A.halfwidth[0]] / (A - B);
      // (The factor of 1/(A-B) insures that the central peak has height 1)

      if (((-filterXY_B.halfwidth[0]<=ix) && (ix<=filterXY_B.halfwidth[0])) &&
          ((-filterXY_B.halfwidth[1]<=iy) && (iy<=filterXY_B.halfwidth[1])))

        filter.aafWeights[iy+halfwidth[1]][ix+halfwidth[0]] -=
          filterXY_B.aafWeights[iy+filterXY_B.halfwidth[1]][ix+filterXY_B.halfwidth[0]] / (A - B);


      //FOR DEBUGGING REMOVE EVENTUALLY
      cerr << "aafWeights["<<iy<<"]["<<ix<<"] = " << filter.aafWeights[iy+halfwidth[1]][ix+halfwidth[0]] << endl;
      //cerr << aafWeights[iy+halfwidth[1]][ix+halfwidth[0]];
      //if (ix == halfwidth[0]) cerr << "\n"; else cerr << " ";

    } // for (int ix=-halfwidth[0]; ix<=halfwidth[0]; ix++) {
  } // for (int iy=-halfwidth[1]; iy<=halfwidth[1]; iy++) {

  cerr << "\n";

  if (pA && pB) {
    *pA = A/(A-B); // Rescale A and B numbers returned to the caller
    *pB = B/(A-B); // (because we divided the array entries by (A-B) earlier)
  }
  return filter;
} //_GenFilterGenDog2D()



template<class RealNum>
// Create a 2-D filter and fill it with a difference of (generalized) Gaussians:
Filter2D<RealNum, int> 
GenFilterGenDog2D(RealNum width_a[2],  //"a" parameter in formula
                  RealNum width_b[2],  //"b" parameter in formula
                  RealNum m_exp,       //"m" parameter in formula
                  RealNum n_exp,       //"n" parameter in formula
                  int halfwidth[2],
                  RealNum *pA=NULL, //optional:report A,B coeffs to user
                  RealNum *pB=NULL) //optional:report A,B coeffs to user
{
  Filter2D<RealNum, int> filterXY_A =
    GenFilterGenGauss2D(width_a,      //"a_x", "a_y" gaussian width parameters
                        m_exp,        //"n" exponent parameter
                        halfwidth);

  Filter2D<RealNum, int> filterXY_B =
    GenFilterGenGauss2D(width_b,      //"b_x", "b_y" gaussian width parameters
                        n_exp,        //"n" exponent parameter
                        halfwidth);

  return _GenFilterGenDog2D(width_a,
                            width_b,
                            m_exp,
                            n_exp,
                            filterXY_A, filterXY_B,
                            pA,
                            pB);
} //GenFilterGenDog2D(...halfwidth...)

template<class RealNum>
// Create a 2-D filter and fill it with a difference of (generalized) Gaussians:
Filter2D<RealNum, int> 
GenFilterGenDog2DThresh(RealNum width_a[2],  //"a" parameter in formula
                        RealNum width_b[2],  //"b" parameter in formula
                        RealNum m_exp,       //"m" parameter in formula
                        RealNum n_exp,       //"n" parameter in formula
                        RealNum window_threshold, //controls filter window width
                        RealNum *pA=NULL, //optional:report A,B coeffs to user
                        RealNum *pB=NULL) //optional:report A,B coeffs to user
{
  Filter2D<RealNum, int> filterXY_A =
    GenFilterGenGauss2DThresh(width_a, //"a_x", "a_y" gaussian width parameters
                              m_exp,   //"m" exponent parameter
                              window_threshold);

  Filter2D<RealNum, int> filterXY_B =
    GenFilterGenGauss2DThresh(width_b, //"b_x", "b_y" gaussian width parameters
                              n_exp,   //"n" exponent parameter
                              window_threshold);

  return _GenFilterGenDog2D(width_a,
                            width_b,
                            m_exp,
                            n_exp,
                            filterXY_A, filterXY_B,
                            pA,
                            pB);
} //GenFilterGenDog2DThresh(...window_threshold...)





template<class RealNum>
RealNum
_ApplyGauss3D(Filter1D<float, int> aFilter[3],
              int const image_size[3], 
              RealNum ***aaafSource,
              RealNum ***aaafDest,
              RealNum ***aaafMask)
              //bool precompute_mask_times_source = true)
{
  assert(aaafSource);
  assert(aaafDest);

  // Optional
  // The "A" 3-D Gaussian coefficient is the product of the
  // 1-D Gaussian coefficients in the X,Y,Z directions.
  // Those coefficients happen to equal the value of the 1-D
  // Gaussian evaluated at its peak, which is stored in the central entry at
  // "halfwidth".  (The 1-D filter arrays have size equal to 2*halfwidth+1)
  RealNum A_coeff = (aFilter[0].afWeights[aFilter[0].halfwidth] *
                     aFilter[1].afWeights[aFilter[1].halfwidth] *
                     aFilter[2].afWeights[aFilter[2].halfwidth]);

  // Create temporary arrays to perform the filter in each direction:
  // (I suppose if I really cared about speed, I could alter Filter1D::Apply() 
  //  function to use pointer arithmatic and allow large strides.
  //  This would also eliminate the need for these temporary arrays.)
  float *aafSource[3];
  float *aafDest[3];
  float *aafMask[3];
  for (int d=0; d < 3; d++) {
    aafDest[d]   = new float [image_size[d]];
    aafSource[d] = new float [image_size[d]];
    aafMask[d]   = NULL;
    if (aaafMask)
      aafMask[d] = new float [image_size[d]];
  }

  // This is a "separable" filter.
  // You can apply the filter sequentially, in the X, Y, Z directions
  // (instead of applying the filter simultaneously in all 3 directions,
  //  which would be much slower)

  // Initially copy aaafSource into aaafDest
  // (We don't want to have to allocate temporary array to 
  //  store the result of each successive filter operation. 
  //  Instead just store the most recent filter operation in aaafDest,
  //  and perform each operation on whatever's currently in aaafDest.)
  for (int iz = 0; iz < image_size[2]; iz++)
    for (int iy = 0; iy < image_size[1]; iy++)
      for (int ix = 0; ix < image_size[0]; ix++)
        aaafDest[iz][iy][ix] = aaafSource[iz][iy][ix];


  int d; //direction where we are applying the filter (x<==>0, y<==>1, z<==>2)


  // Apply the filter in the Z direction (d=2):
  d = 2;
  cerr << "  progress: Applying Z filter. Processing Y plane#" << endl;
  for (int iy = 0; iy < image_size[1]; iy++) {
    cerr << "  " << iy+1 << " / " << image_size[1] << "\n";
    for (int ix = 0; ix < image_size[0]; ix++) {
      // copy the data we need to the temporary array
      for (int iz = 0; iz < image_size[2]; iz++) {
        aafSource[d][iz] = aaafDest[iz][iy][ix];  // copy from aaafDest
        if (aaafMask)
          aafMask[d][iz] = aaafMask[iz][iy][ix];
      }
      // apply the filter to the 1-D temporary array
      aFilter[d].Apply(image_size[d],
                       aafSource[d],
                       aafDest[d],
                       aafMask[d]);
                       //precompute_mask_times_source);
      for (int iz = 0; iz < image_size[d]; iz++)
        aaafDest[iz][iy][ix] = aafDest[d][iz];  // copy back into aaafDest
    } //for (int ix = 0; ix < image_size[0]; ix++)
  } //for (int iy = 0; iy < image_size[1]; iy++)


  // Apply the filter in the Y direction:
  d=1;
  cerr << "  progress: Applying Y filter. Processing Z plane#" << endl;
  for (int iz = 0; iz < image_size[2]; iz++) {
    cerr << "  " << iz+1 << " / " << image_size[2] << "\n";
    for (int ix = 0; ix < image_size[0]; ix++) {
      // copy the data we need to the temporary array
      for (int iy = 0; iy < image_size[1]; iy++) {
        aafSource[d][iy] = aaafDest[iz][iy][ix]; //data from previous aaafDest
        if (aaafMask)
          aafMask[d][iy] = aaafMask[iz][iy][ix];
      }
      // apply the filter to the 1-D temporary array
      aFilter[d].Apply(image_size[d],
                       aafSource[d],
                       aafDest[d],
                       aafMask[d]);
                       //precompute_mask_times_source);
      for (int iy = 0; iy < image_size[d]; iy++)
        aaafDest[iz][iy][ix] = aafDest[d][iy];  // copy back into aaafDest
    } //for (int ix = 0; ix < image_size[0]; ix++) {
  } //for (int iz = 0; iz < image_size[2]; iz++)

  // Apply the filter in the X direction:
  d=0;
  cerr << "  progress: Applying X filter. Processing Z plane#" << endl;
  for (int iz = 0; iz < image_size[2]; iz++) {
    cerr << "  " << iz+1 << " / " << image_size[2] << "\n";
    for (int iy = 0; iy < image_size[1]; iy++) {
      // copy the data we need to the temporary array
      for (int ix = 0; ix < image_size[0]; ix++) {
        aafSource[d][ix] = aaafDest[iz][iy][ix]; //data from previous aaafDest
        if (aaafMask)
          aafMask[d][ix] = aaafMask[iz][iy][ix];
      }
      // apply the filter to the 1-D temporary array
      aFilter[d].Apply(image_size[d],
                       aafSource[d],
                       aafDest[d],
                       aafMask[d]);
                       //precompute_mask_times_source);
      for (int ix = 0; ix < image_size[d]; ix++)
        aaafDest[iz][iy][ix] = aafDest[d][ix];  // copy back into aaafDest
    } //for (int iy = 0; iy < image_size[1]; iy++)
  } //for (int iz = 0; iz < image_size[2]; iz++)

  // delete the temporary arrays
  for (int d=0; d<3; d++) {
    delete [] aafSource[d];
    delete [] aafDest[d];
    if (aafMask[d])
      delete [] aafMask[d];
  }

  return A_coeff;
} //_ApplyGauss3D(aFilter)




template<class RealNum>
RealNum
ApplyGauss3D(RealNum const sigma[3],
             int const window_halfwidth[3],
             int const image_size[3], 
             RealNum ***aaafSource,
             RealNum ***aaafDest,
             RealNum ***aaafMask)
             //bool precompute_mask_times_source = true)
{
  assert(aaafSource);
  assert(aaafDest);
  //assert(aaafMask);
  //Allocate filters in all 3 directions.  (Later apply them sequentially.)
  Filter1D<float, int> aFilter[3];
  for (int d=0; d < 3; d++)
    aFilter[d] = GenFilterGauss1D(sigma[d], window_halfwidth[d]);

  return _ApplyGauss3D(aFilter,
                       image_size, 
                       aaafSource,
                       aaafDest,
                       aaafMask);
                       //precompute_mask_times_source)

} //ApplyGauss3D(sigma, window_halfwidth, ...)


template<class RealNum>
RealNum
ApplyGauss3D(RealNum const sigma[3],
             RealNum window_threshold,
             int const image_size[3], 
             RealNum ***aaafSource,
             RealNum ***aaafDest,
             RealNum ***aaafMask)
             //bool precompute_mask_times_source = true)
{
  assert(aaafSource);
  assert(aaafDest);
  //assert(aaafMask);
  //Allocate filters in all 3 directions.  (Later apply them sequentially.)
  Filter1D<float, int> aFilter[3];
  for (int d=0; d < 3; d++)
    aFilter[d] = GenFilterGauss1DThresh(sigma[d], window_threshold);

  return _ApplyGauss3D(aFilter,
                       image_size, 
                       aaafSource,
                       aaafDest,
                       aaafMask);
                       //precompute_mask_times_source)
} //ApplyGauss3D(sigma, window_threshold, ...)




int main(int argc, char **argv) {
  try {
    Settings settings; // parse the command-line argument list from the shell
    settings.ParseArgs(argc, argv);

    // Read the input tomogram
    cerr << "Reading tomogram \""<<settings.in_file_name<<"\"" << endl;
    MrcSimple tomo_in;
    tomo_in.Read(settings.in_file_name, false);
    // (Note: You can also use "tomo_in.Read(cin);" or "cin >> tomo;")
    tomo_in.PrintStats(cerr);      //Optional (display the tomogram size & format)

    // ---- mask ----

    // Optional: if there is a "mask", read that too
    MrcSimple mask;
    if (settings.mask_file_name != "") {
      cerr << "Reading mask \""<<settings.mask_file_name<<"\"" << endl;
      mask.Read(settings.mask_file_name, false);
      if ((mask.mrc_header.nvoxels[0] != tomo_in.mrc_header.nvoxels[0]) ||
          (mask.mrc_header.nvoxels[1] != tomo_in.mrc_header.nvoxels[1]) ||
          (mask.mrc_header.nvoxels[2] != tomo_in.mrc_header.nvoxels[2]))
        throw InputErr("Error: The size of the mask does not match the size of the tomogram.\n");
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
      tomo_in.Rescale01(mask.aaafDensity);

    // ---- make an array that will store the new tomogram we will create ----

    cerr << "allocating space for new tomogram..." << endl;
    MrcSimple tomo_out = tomo_in; //this will take care of allocating the array

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
      voxel_width[0] = tomo_in.mrc_header.cellA[0]/tomo_in.mrc_header.nvoxels[0];
      voxel_width[1] = tomo_in.mrc_header.cellA[1]/tomo_in.mrc_header.nvoxels[1];
      voxel_width[2] = tomo_in.mrc_header.cellA[2]/tomo_in.mrc_header.nvoxels[2];
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
      throw InputErr("Error in tomogram header: Invalid voxel width(s).\n"
                     "Use the -w argument to specify the voxel width in nm.");

    if (abs((voxel_width[0] - voxel_width[1])
            /
            (0.5*(voxel_width[0] + voxel_width[1]))) > 0.0001)
      throw InputErr("Error in tomogram header: Unequal voxel widths in the x and y directions.\n"
                     "Use the -w argument to specify the voxel width in nm.");
    for (int d=0; d<3; d++) {
      settings.width_a[d] /= voxel_width[d];
      settings.width_b[d] /= voxel_width[d];
      //settings.window_halfwidth[d] = floor(settings.window_ratio * 
      //                                    MAX(settings.width_a[d],
      //                                        settings.width_b[d]));
    }

    //cerr << "applying filter (window size in voxels: "
    //     << 1 + 2*settings.window_halfwidth[0] << ","
    //     << 1 + 2*settings.window_halfwidth[1] << ","
    //     << 1 + 2*settings.window_halfwidth[2] << ")"
    //     << " ..." << endl;

    if (settings.filter_type == Settings::NONE) {
      cerr << "filter_type = Intensity Map <No convolution filter specified>\n";
      // Not needed:
      //for (int iz = 0; iz < size[2]; iz++)
      //  for (int iy = 0; iy < size[1]; iy++)
      //    for (int ix = 0; ix < size[0]; ix++)
      //      tomo_out.aaafDensity[iz][iy][ix]=tomo_in.aaafDensity[iz][iy][ix];
      // (We have copied the contents from tomo_in into tomo_out already.)
    } 

    else if (settings.filter_type == Settings::GAUSS) {
      cerr << "filter_type = Gaussian\n";

      float A; // let the user know what A coefficient was used

      int window_halfwidth[3] = {-1, -1, -1};
      if (settings.window_ratio > 0)
        for (int d=0; d < 3; d++)
          window_halfwidth[d] = floor(settings.width_a[d] *
                                      settings.window_ratio);

      // Convolve the original source with the 1st Gaussian
      if (window_halfwidth[0] > 0)
        A = ApplyGauss3D(settings.width_a,
                         window_halfwidth,
                         tomo_in.mrc_header.nvoxels,
                         tomo_in.aaafDensity,
                         tomo_out.aaafDensity,
                         mask.aaafDensity);
                         //true);
      else if (settings.window_threshold > 0)
        A = ApplyGauss3D(settings.width_a,
                         settings.window_threshold,
                         tomo_in.mrc_header.nvoxels,
                         tomo_in.aaafDensity,
                         tomo_out.aaafDensity,
                         mask.aaafDensity);
                         //true);
      else
        assert(false);

      cerr << " Filter Used:\n"
        " h(x,y,z)   = A*exp(-0.5*((x/a_x)^2 + (y/a_y)^2 + (z/a_z)^))\n"
        " ... where  A,  a_x (in voxels)  equal:\n"
           << " " << A << " " 
           << " " << settings.width_a[0] << "\n"
        " ...   and  A,  a_y (in voxels)  equal:\n"
           << " " << A << " "
           << " " << settings.width_a[1] << "\n"
        " ...   and  A,  a_z (in voxels)  equal:\n"
           << " " << A << " "
           << " " << settings.width_a[2] << "\n"
          " You can plot this function in the X,Y, or Z directions using:\n"
          " draw_filter_1D.py -gauss  A  a\n";
        
    } // if (settings.filter_type == Settings::GAUSS)


    else if (settings.filter_type == Settings::GGAUSS) {
      throw InputErr("Error: 3-D Generalized-Gaussians filter\n"
                     "       has not been implemented yet.\n"
                     "       (However implementing should be trivial to do.)\n"
                     "For now, you must avoid using the -exponent argument.  Use ordinary Gaussians\n");

      cerr << " Filter Used:\n"
        " h(x,y,z)   = A*exp(-((x/a_x)^2 + (y/a_y)^2 + (z/a_z)^2)^(m/2))\n"
        " ... where  A,  a_x(voxels),  m  equal:\n"
           << " " << A << " " 
           << " " << settings.width_a[0]
           << " " << settings.exp_m << "\n"
        " ...   and  A,  a_y(voxels),  m  equal:\n"
           << " " << A << " "
           << " " << settings.width_a[1]
           << " " << settings.exp_m << "\n"
        " ...   and  A,  a_z(voxels),  m  equal:\n"
           << " " << A << " "
           << " " << settings.width_a[2]
           << " " << settings.exp_m << "\n"
        " You can plot this function in the X,Y, or Z directions using:\n"
        " draw_filter_1D.py -ggauss  A  a  m\n";
    }


    else if (settings.filter_type == Settings::DOG) {
      cerr << "filter_type = Difference of Gaussians (DOG)\n";

      float ***aaafTemp; //temporary array to store partially processed tomogram
      float *afTemp;     //temporary array to store partially processed tomogram
      int   *size = tomo_in.mrc_header.nvoxels;
      Alloc3D(tomo_in.mrc_header.nvoxels,
              &afTemp,
              &aaafTemp);

      float A, B;        // let the user know what A B coefficients were used

      int window_halfwidth[3] = {-1, -1, -1};
      if (settings.window_ratio > 0)
        for (int d=0; d < 3; d++)
          window_halfwidth[d] = floor(settings.window_ratio *
                                      MAX(settings.width_a[d],
                                          settings.width_b[d]));
                                     

      // Convolve the original source with the 1st Gaussian
      if (window_halfwidth[0] > 0)
        A = ApplyGauss3D(settings.width_a,
                         window_halfwidth,
                         tomo_in.mrc_header.nvoxels,
                         tomo_in.aaafDensity,
                         tomo_out.aaafDensity,
                         mask.aaafDensity);
                         //true);
      else if (settings.window_threshold > 0)
        A = ApplyGauss3D(settings.width_a,
                         settings.window_threshold,
                         tomo_in.mrc_header.nvoxels,
                         tomo_in.aaafDensity,
                         tomo_out.aaafDensity,
                         mask.aaafDensity);
                         //true);
      else
        assert(false);

      // Convolve the original source with the 2nd Gaussian
      if (window_halfwidth[0] > 0)
        B = ApplyGauss3D(settings.width_b,
                         window_halfwidth,
                         tomo_in.mrc_header.nvoxels,
                         tomo_in.aaafDensity,
                         aaafTemp,
                         mask.aaafDensity);
                         //true);
      else if (settings.window_threshold > 0)
        B = ApplyGauss3D(settings.width_b,
                         settings.window_threshold,
                         tomo_in.mrc_header.nvoxels,
                         tomo_in.aaafDensity,
                         aaafTemp,
                         mask.aaafDensity);
                         //true);
      else
        assert(false);

      for (int iz = 0; iz < size[2]; iz++) {
        for (int iy = 0; iy < size[1]; iy++) {
          for (int ix = 0; ix < size[0]; ix++) {
            // Subtract the second convolved signal from the first
            tomo_out.aaafDensity[iz][iy][ix] -= aaafTemp[iz][iy][ix];
            // Optional: Scale by 1/(A-B) so that the central peak has height 1
            tomo_out.aaafDensity[iz][iy][ix] *= (1.0 / (A-B));
          }
        }
      }

      float Aeff = A/(A-B); // rescale so that central peak has height 1
      float Beff = B/(A-B); // rescale so that central peak has height 1
      cerr << " Filter Used:\n"
        " h(x,y,z)   = h_a(x,y,z) - h_b(x,y,z)\n"
        " h_a(x,y,z) = A*exp(-0.5*((x/a_x)^2 + (y/a_y)^2 + (z/a_z)^2))\n"
        " h_b(x,y,z) = B*exp(-0.5*((x/b_x)^2 + (y/b_y)^2 + (z/b_z)^2))\n"
        " ... where  A,  B,  a_x,  b_x  (in voxels) equal:\n"
           << " " << Aeff << " " << Beff
           << " " << settings.width_a[0]
           << " " << settings.width_b[0] << "\n"
        " ...   and A,  B,  a_y,  b_y  (in voxels) equal:\n"
           << " " << Aeff << " " << Beff
           << " " << settings.width_a[1]
           << " " << settings.width_b[1] << "\n"
        " ...   and A,  B,  a_z,  b_z  (in voxels) equal:\n"
           << " " << Aeff << " " << Beff
           << " " << settings.width_a[2]
           << " " << settings.width_b[2] << "\n"
          " You can plot this function in the X,Y, or Z directions using:\n"
          " draw_filter_1D.py  -dog  A  B  a  b\n";                            

      // Deallocate the temporary array
      Dealloc3D(tomo_in.mrc_header.nvoxels,
                &afTemp,
                &aaafTemp);

    } //if (settings.filter_type == Settings::DOG)

    else if (settings.filter_type == Settings::DOGG) {
      cerr << "filter_type = Difference-of-Generalized-Gaussians (DOGG)\n";

      throw InputErr("Error: The 3-D Difference-of-Generalized-Gaussians filter\n"
                     "       has not been implemented yet.\n"
                     "      (However implementing should be trivial to do.\n"
                     "       Edit the code to define a functions \"GenFilterGenDog3D()\" and\n"
                     "       GenFilterGenGauss3D This should be easy. They will be nearly identical\n"
                     "       to the functions \"GenFilterGenDog2D()\" and \"GenFilterGenGauss2D()\")\n"
                     "For now, you can use an ordinary 3-D DOG filter with default exponents m=n=2.\n");

      //Filter3D<float, int> filterXY = 
      //GenFilterGenDog3D(settings.width_a,  //"a" parameter in formula
      //                  settings.width_b,  //"b" parameter in formula
      //                  settings.m_exp,    //"m" parameter in formula
      //                  settings.n_exp,    //"n" parameter in formula
      //                  settings.window_threshold,
      //                  &A,
      //                  &B);
      cerr << " Filter Used:\n"
        " h(x,y,z)   = h_a(x,y,z) - h_b(x,y,z)\n"
        " h_a(x,y,z) = A*exp(-((x/a_x)^2 + (y/a_y)^2 + (z/a_z)^2)^(m/2))\n"
        " h_b(x,y,z) = B*exp(-((x/b_x)^2 + (y/b_y)^2 + (z/b_z)^2)^(n/2))\n"
        " ... where, in the X direction,  A,  B,  a_x,  b_x(in voxels),  m,  n  equal:\n"
           << " " << A << " " << B
           << " " << settings.width_a[0]
           << " " << settings.width_b[0]
           << " " << settings.m_exp
           << " " << settings.n_exp << "\n"
        " ...    and in the Y direction,  A,  B,  a_y,  b_y(in voxels),  m,  n  equal:\n"
           << " " << A << " " << B
           << " " << settings.width_a[1]
           << " " << settings.width_b[1]
           << " " << settings.m_exp
           << " " << settings.n_exp << "\n"
        " ...    and in the Z direction,  A,  B,  a_z,  b_z(in voxels),  m,  n  equal:\n"
           << " " << A << " " << B
           << " " << settings.width_a[2]
           << " " << settings.width_b[2]
           << " " << settings.m_exp
           << " " << settings.n_exp << "\n"
          " You can plot these functions using:\n"
          " draw_filter_1D.py  -dogg  A  B  a  b  m  n\n";
    } //else if (settings.filter_type == Settings::DOGG)


    else if (settings.filter_type == Settings::DOGGXY) {
      cerr << "filter_type = Difference-of-Generalized-Gaussians in the XY plane\n";
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
      //This reduces the computation by a factor of O(filter.halfwidth[2]))
      //(A 1-D convolution followed by a 2-D convolution is much faster per 
      // voxel than a full 3-D convolution.)

      // First, generate the filter in the Z direction:

      int window_halfwidthZ = -1;
      if (settings.window_ratio > 0.0)
        window_halfwidthZ = floor(settings.width_a[2] *
                                  settings.window_ratio);

      Filter1D<float, int> filterZ;
      if (window_halfwidthZ > 0)
        filterZ = GenFilterGauss1D(settings.width_a[2],
                                   window_halfwidthZ);
      else if (settings.window_threshold > 0.0) {
        window_halfwidthZ = floor(settings.width_a[2] *
                                  sqrt(-2*log(settings.window_threshold)));
        filterZ = GenFilterGauss1D(settings.width_a[2],
                                   window_halfwidthZ);
      }
      float C;       // let the user know what C coefficient was used
      C = filterZ.afWeights[window_halfwidthZ]; //(C=peak height located at the
                                                //   middle of the array)


      // Then generate the filter in the XY directions

      Filter2D<float, int> filterXY;
      float A, B;       // let the user know what A B coefficients were used

      if (settings.window_ratio > 0.0) {
        // If the user did not specify the width of the filters explicitly,
        // then, determine the filter window size (filter width) from the
        // "window_ratio" parameters which are distances expressed
        // as multiples of the widths of the gaussians (a and b parameters)
        int window_halfwidth[2];
        for (int d=0; d < 2; d++)
          window_halfwidth[d] = floor(settings.window_ratio *
                                      MAX(settings.width_a[d],
                                          settings.width_b[d]));;
                                     
        filterXY = GenFilterGenDog2D(settings.width_a,//"a" parameter in formula
                                     settings.width_b,//"b" parameter in formula
                                     settings.m_exp,  //"m" parameter in formula
                                     settings.n_exp,  //"n" parameter in formula
                                     window_halfwidth,
                                     &A,
                                     &B);
      } //else if (settings.window_ratio > 0.0)
      else if (settings.window_threshold > 0)
        filterXY = GenFilterGenDog2DThresh(settings.width_a,//"a" parameter
                                     settings.width_b,//"b" parameter in formula
                                     settings.m_exp,  //"m" parameter in formula
                                     settings.n_exp,  //"n" parameter in formula
                                     settings.window_threshold,
                                     &A,
                                     &B);
      else
        assert(false);      

      // Precompute the effect of the filter in the Z direction.

      // Create temporary 1-D arrays to perform the filter in the Z-direction:
      float *afDensityZorig = new float [tomo_in.mrc_header.nvoxels[2]];
      float *afDensityZnew  = new float [tomo_in.mrc_header.nvoxels[2]];
      float *afMask         = NULL;
      if (mask.aaafDensity)
        afMask       = new float [tomo_in.mrc_header.nvoxels[2]];

      // Then apply the filter in the Z direction
      // (and store the filtered 3-D image in the original tomogram array)
      for (int ix = 0; ix < tomo_in.mrc_header.nvoxels[0]; ix++) {
        for (int iy = 0; iy < tomo_in.mrc_header.nvoxels[1]; iy++) {
          for (int iz = 0; iz < tomo_in.mrc_header.nvoxels[2]; iz++) {
            afDensityZorig[iz] = tomo_in.aaafDensity[iz][iy][ix];
            if (afMask)
              afMask[iz] = mask.aaafDensity[iz][iy][ix];
          }
          filterZ.Apply(tomo_in.mrc_header.nvoxels[2],
                        afDensityZorig,
                        afDensityZnew,
                        afMask);
                        //true);
          //It would be wasteful to allocate a temporary array to store this
          //Instead store the result of the convolution in the original array:
          for (int iz = 0; iz < tomo_in.mrc_header.nvoxels[2]; iz++)
            tomo_in.aaafDensity[iz][iy][ix] = afDensityZnew[iz];
        } //for (int iy = 0; iy < tomo_in.mrc_header.nvoxels[1]; iy++) {
      } // for (int ix = 0; ix < tomo_in.mrc_header.nvoxels[0]; ix++) {
      delete [] afDensityZorig;
      delete [] afDensityZnew;
      if (afMask)
        delete [] afMask;

      // Now apply the filter in the X and Y directions:
      cerr << "  progress: processing plane#" << endl;
      for (int iz = 0; iz < tomo_in.mrc_header.nvoxels[2]; iz++) {
        cerr << "  " << iz+1 << " / " << tomo_in.mrc_header.nvoxels[2] << "\n";
        float **aafMaskXY = NULL;
        if (mask.aaafDensity)
          aafMaskXY = mask.aaafDensity[iz];
        filterXY.Apply(tomo_in.mrc_header.nvoxels,
                       tomo_in.aaafDensity[iz],
                       tomo_out.aaafDensity[iz],
                       tomo_out.aaafDensity[iz]);
                       //aafMaskXY);
                       //true);
      }

      cerr << " Filter Used:\n"
        " h(x,y,z) = (h_a(x,y) - h_b(x,y)) * C * exp(-0.5*(z/s)^2)\n"
        " h_a(x,y) = A*exp(-((x/a_x)^2 + (y/a_y)^2)^(m/2))\n"
        " h_b(x,y) = B*exp(-((x/b_x)^2 + (y/b_y)^2)^(n/2))\n"
        " ... where, in the X direction,  A,  B,  a_x,  b_x(in voxels),  m,  n  equal:\n"
           << " " << A << " " << B
           << " " << settings.width_a[0]
           << " " << settings.width_b[0]
           << " " << settings.m_exp
           << " " << settings.n_exp << "\n"
        " ...    and in the Y direction,  A,  B,  a_y,  b_y(in voxels),  m,  n  equal:\n"
           << " " << A << " " << B
           << " " << settings.width_a[1]
           << " " << settings.width_b[1]
           << " " << settings.m_exp
           << " " << settings.n_exp << "\n"
        " You can plot these functions using:\n"
        " draw_filter_1D.py  -dogg  A  B  a  b  m  n\n"
        "\n"
        " ...    and in the Z direction,  C, s(in voxels)  equals:\n"
           << " " << C <<
           << " " << settings.width_a[2] << "\n"
          " You can plot this function using:\n"
          " draw_filter_1D.py  -gauss C s\n";

  
    } //else if (settings.filter_type = Settings::DOGGXY)





    // --- Exchange light voxels for dark voxels ? ---

    if (settings.invert_output)
      tomo_out.Invert(mask.aaafDensity);


    // --- Rescale so that the lowest, highest voxels have density 0 and 1? ---

    if (settings.in_rescale01)
      tomo_out.Rescale01(mask.aaafDensity);





    // ----- thresholding and masking: -----
    
    if (settings.use_thresholds) {

      cerr << "thresholding..." << endl;

      for (int iz=0; iz<tomo_out.mrc_header.nvoxels[2]; iz++) {
        for (int iy=0; iy<tomo_out.mrc_header.nvoxels[1]; iy++) {
          for (int ix=0; ix<tomo_out.mrc_header.nvoxels[0]; ix++) {
            if (! settings.use_dual_thresholds)
              tomo_out.aaafDensity[iz][iy][ix] =
                Threshold2(tomo_out.aaafDensity[iz][iy][ix],
                           settings.out_threshold_01_a,
                           settings.out_threshold_01_b);
            else
              tomo_out.aaafDensity[iz][iy][ix] =
                Threshold4(tomo_out.aaafDensity[iz][iy][ix],
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
              tomo_out.aaafDensity[iz][iy][ix] = settings.mask_out;


    // Write the file
    if (settings.out_file_name != "") {
      cerr << "writing tomogram (in float mode)" << endl;
      tomo_out.Write(settings.out_file_name);
      //(You can also use "tomo_out.Write(cout);" or "cout<<tomo_out;")
    }
  }

  catch (InputErr& e) {
    cerr << "\n" << e.what() << endl;
    exit(1);
  }
}







