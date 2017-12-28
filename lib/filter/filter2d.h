#ifndef _FILTER2D_H
#define _FILTER2D_H

#include <ostream>
using namespace std;





template<class RealNum, class Integer>

class Filter2D {
public:
  RealNum *af;
  RealNum **aafWeights;
  Integer size[2];



  // Apply the filter to tomographic data in the "aafSource" array
  // Save the results in the "aafDest" array.  (A "mask" is optional.)
  // All arrays are 2D and assumed to be the same size, given by size_source[].
  void Apply(Integer size_source[2],
             RealNum **aafSource,
             RealNum **aafDest,
             RealNum **aafMask = NULL,
             bool precompute_mask_times_source = false) const
  {

    // Apply the filter to the original tomogram data. (Store in aafSource)
    //        ___
    //        \
    // g(i) = /__  h(j) * f(i-j)
    //         j
    //
    // where: f(i) is the original density of the tomogram at position ix,iy
    //        h(j) is the filter (smoothing function)
    //       
    //  Note on mask functions:
    //     When summing over j, we ignore contributions from voxels where
    //     mask(i-j) is zero.  We don't count them in the average.
    //     Because h(j) is not necessarily normalized, g(i) is later divided by
    //     the area under the curve h(j)*mask(i-j)    (as a function of j)
    //     

    if (aafMask && precompute_mask_times_source) {
      // The mask should be 1 everywhere we want to consider, and 0 elsewhere.
      // Multiplying the density in the tomogram by the mask removes some of 
      // the voxels from consideration later on when we do the filtering.
      // (Later, we will adjust the weight of the average we compute when we
      //  apply the filter in order to account for the voxels we deleted now.)
      for (int iy=0; iy<size_source[1]; iy++)
        for (int ix=0; ix<size_source[0]; ix++)
          aafSource[iy][ix] *= aafMask[iy][ix];
    }

    for (Integer iy=0; iy<size_source[1]; iy++) {

      for (Integer ix=0; ix<size_source[0]; ix++) {

        if ((aafMask) && (aafMask[iy][ix] == 0.0))
          continue;
          
        RealNum g = 0.0;
        RealNum denominator = 0.0;

        for (Integer jy=-size[1]; jy<=size[1]; jy++) {

          if ((iy-jy < 0) || (size_source[1] <= iy-jy))
            continue;

          for (Integer jx=-size[0]; jx<=size[0]; jx++) {

            if ((ix-jx < 0) || (size_source[0] <= ix-jx))
              continue;

            RealNum delta_g = 
              aafWeights[jy+size[1]][jx+size[0]] * aafSource[iy-jy][ix-jx];

            if (! precompute_mask_times_source)
              delta_g *= aafMask[iy-jy][ix-jx];

            g += delta_g;
              // Note: We previously applied the mask by multiplying
              //          aaDensity[iy][ix]
              //         by aafMask[iy][ix]   (if present)
            if (aafMask)
              denominator +=
                aafMask[iy-jy][ix-jx] * aafWeights[jy+size[1]][jx+size[0]];
                                               
            else
              denominator += aafWeights[jy+size[1]][jx+size[0]];
                                          
            // Note: If there were no mask, and if the filter is normalized
            // then denominator=1 always, and we could skip the line above.
          }
        }


        if (denominator > 0.0)
          g /= denominator;
        else
          //Otherwise, this position lies outside the mask region.
          g = 0.0;

        aafDest[iy][ix] = g;
      }
    }
  } //Filter2D::Apply()


  Filter2D(Integer size_half[2]) {
    Integer table_size[2];
    for(Integer d=0; d < 2; d++) {
      size[d] = size_half[d];
      table_size[d] = 1 + 2*size_half[d];
    }
    Alloc2D(table_size, &af, &aafWeights);
  }
  ~Filter2D() {
    Integer table_size[2];
    for(Integer d=0; d < 2; d++)
      table_size[d] = 1 + 2*size[d];
    Dealloc2D(table_size, &af, &aafWeights);
  }
  void Normalize() {
    // Make sure the sum of the filter weights is 1
    RealNum total = 0.0;
    for (Integer iy=-size[1]; iy<=size[1]; iy++)
      for (Integer ix=-size[0]; ix<=size[0]; ix++)
        total += aafWeights[iy+size[1]][ix+size[0]];
    for (Integer iy=-size[1]; iy<=size[1]; iy++)
      for (Integer ix=-size[0]; ix<=size[0]; ix++)
        aafWeights[iy+size[1]][ix+size[0]] /= total;
  }
}; // class Filter2D



#endif //#ifndef _FILTER2D_H
