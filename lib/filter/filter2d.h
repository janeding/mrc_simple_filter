#ifndef _FILTER2D_H
#define _FILTER2D_H

#include <cstring>
#include <ostream>
using namespace std;





template<class RealNum, class Integer>

class Filter2D {
public:
  RealNum *af;
  RealNum **aafWeights;
  Integer halfwidth[2]; //num pixels from filter center to edge in x,y directions



  // Apply the filter to tomographic data in the "aafSource" array
  // Save the results in the "aafDest" array.  (A "mask" is optional.)
  // All arrays are 2D and assumed to be the same size
  void Apply(Integer const size_source[2],
             RealNum **aafSource,
             RealNum **aafDest,
             RealNum **aafMask = NULL) const
             //bool precompute_mask_times_source = true) const
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

    // The mask should be 1 everywhere we want to consider, and 0 elsewhere.
    // Multiplying the density in the tomogram by the mask removes some of 
    // the voxels from consideration later on when we do the filtering.
    // (Later, we will adjust the weight of the average we compute when we
    //  apply the filter in order to account for the voxels we deleted now.)
    // Precomupting the mask product is faster but modifies the source image.

    //if (aafMask && precompute_mask_times_source)
    if (aafMask) 
      for (int iy=0; iy<size_source[1]; iy++)
        for (int ix=0; ix<size_source[0]; ix++)
          aafSource[iy][ix] *= aafMask[iy][ix];

    for (Integer iy=0; iy<size_source[1]; iy++) {

      for (Integer ix=0; ix<size_source[0]; ix++) {

        if ((aafMask) && (aafMask[iy][ix] == 0.0))
          continue;
          
        RealNum g = 0.0;
        RealNum denominator = 0.0;

        for (Integer jy=-halfwidth[1]; jy<=halfwidth[1]; jy++) {

          if ((iy-jy < 0) || (size_source[1] <= iy-jy))
            continue;

          for (Integer jx=-halfwidth[0]; jx<=halfwidth[0]; jx++) {

            if ((ix-jx < 0) || (size_source[0] <= ix-jx))
              continue;

            RealNum delta_g = 
              aafWeights[jy+halfwidth[1]][jx+halfwidth[0]] * aafSource[iy-jy][ix-jx];

            //if (! precompute_mask_times_source)
            //  delta_g *= aafMask[iy-jy][ix-jx];

            g += delta_g;
              // Note: We previously applied the mask by multiplying
              //          aaDensity[iy][ix]
              //         by aafMask[iy][ix]   (if present)
            if (aafMask)
              denominator +=
                aafMask[iy-jy][ix-jx] * aafWeights[jy+halfwidth[1]][jx+halfwidth[0]];
                                               
            else
              denominator += aafWeights[jy+halfwidth[1]][jx+halfwidth[0]];
                                          
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


  void Alloc(Integer const set_halfwidth[2]) {
    Integer array_size[2];
    for(Integer d=0; d < 2; d++) {
      halfwidth[d] = set_halfwidth[d];
      array_size[d] = 1 + 2*halfwidth[d];
    }
    Alloc2D(array_size, &af, &aafWeights);
  }

  void Dealloc() {
    Integer array_size[2];
    for(Integer d=0; d < 2; d++) {
      array_size[d] = 1 + 2*halfwidth[d];
      halfwidth[d] = -1;
    }
    Dealloc2D(array_size, &af, &aafWeights);
  }


  void Resize(Integer const set_halfwidth[2]) {
    Integer array_size[2];
    if (af && aafWeights) {
      for(Integer d=0; d < 2; d++)
        array_size[d] = 1 + 2*halfwidth[d];
      Dealloc2D(array_size, &af, &aafWeights);
    }
    for(Integer d=0; d < 2; d++) {
      halfwidth[d] = set_halfwidth[d];
      array_size[d] = 1 + 2*halfwidth[d];
    }
    Alloc2D(array_size, &af, &aafWeights);
  }

  Filter2D(Integer const set_halfwidth[2]) {
    af = NULL;
    aafWeights = NULL;
    Resize(set_halfwidth);
  }

  Filter2D() {
    halfwidth[0] = -1;
    halfwidth[1] = -1;
    af = NULL;
    aafWeights = NULL;
  }

  ~Filter2D() {
    Dealloc();
  }

  void Normalize() {
    // Make sure the sum of the filter weights is 1
    RealNum total = 0.0;
    for (Integer iy=-halfwidth[1]; iy<=halfwidth[1]; iy++)
      for (Integer ix=-halfwidth[0]; ix<=halfwidth[0]; ix++)
        total += aafWeights[iy+halfwidth[1]][ix+halfwidth[0]];
    for (Integer iy=-halfwidth[1]; iy<=halfwidth[1]; iy++)
      for (Integer ix=-halfwidth[0]; ix<=halfwidth[0]; ix++)
        aafWeights[iy+halfwidth[1]][ix+halfwidth[0]] /= total;
  }

  inline Filter2D<RealNum, Integer>&
    operator = (const Filter2D<RealNum, Integer>& source) {
    Resize(source.halfwidth); // allocates and initializes af and aaafWeights
    //for(Int iy=-halfwidth[1]; iy<=halfwidth[1]; iy++)
    //  for(Int ix=-halfwidth[0]; ix<=halfwidth[0]; ix++)
    //    aafWeights[iy][ix] = source.aafWeights[iy][ix];
    // Use memcpy() instead:
    memcpy(af,
           source.af,
           ((1+2*halfwidth[0]) * (1+2*halfwidth[1]))
           *sizeof(RealNum));
  } // operator = ()

}; // class Filter2D



#endif //#ifndef _FILTER2D_H
