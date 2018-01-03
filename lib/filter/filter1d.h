#ifndef _FILTER1D_H
#define _FILTER1D_H

using namespace std;




template<class RealNum, class Integer>

class Filter1D {
public:
  RealNum *afWeights;
  Integer size;

  // Apply the filter to tomographic data in the "afSource" array
  // Save the results in the "afDest" array.  (A "mask" is optional.)
  // All arrays are 1D and assumed to be the same size, given by size_source.
  void Apply(Integer size_source,
             RealNum **afSource,
             RealNum **afDest,
             RealNum **afMask = NULL,
             bool precompute_mask_times_source = false) const
  {

    // Apply the filter to the original tomogram data. (Store in afSource)
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

    if (afMask && precompute_mask_times_source) {
      // The mask should be 1 everywhere we want to consider, and 0 elsewhere.
      // Multiplying the density in the tomogram by the mask removes some of 
      // the voxels from consideration later on when we do the filtering.
      // (Later, we will adjust the weight of the average we compute when we
      //  apply the filter in order to account for the voxels we deleted now.)
      for (int ix=0; ix<size_source; ix++)
        afSource[ix] *= afMask[ix];
    }

    for (Integer ix=0; ix<size_source; ix++) {

      if ((afMask) && (afMask[ix] == 0.0))
        continue;
          
      RealNum g = 0.0;
      RealNum denominator = 0.0;

      for (Integer jx=-size; jx<=size; jx++) {

        if ((ix-jx < 0) || (size_source <= ix-jx))
          continue;

        RealNum delta_g = afWeights[jx+size] * afSource[ix-jx];

        if (! precompute_mask_times_source)
          delta_g *= afMask[ix-jx];

        g += delta_g;
          // Note: We previously applied the mask by multiplying
          //          aDensity[ix]
          //         by afMask[ix]   (if present)
        if (afMask)
          denominator += afMask[ix-jx] * afWeights[jx+size];
        else
          denominator += afWeights[jx+size];
                                          
        // Note: If there were no mask, and if the filter is normalized
        // then denominator=1 always, and we could skip the line above.
      }

      if (denominator > 0.0)
        g /= denominator;
      else
        //Otherwise, this position lies outside the mask region.
        g = 0.0;

      afDest[ix] = g;
    }
  } //Filter1D::Apply()

  void Resize(Integer size_half) {
    if (afWeights)
      delete [] afWeights;
    size = size_half;
    Integer table_size = 1 + 2*size_half;
    afWeights = new RealNum [table_size];
  }

  Filter1D() {
    size = -1;
    afWeights = NULL;
  }

  Filter1D(Integer size_half) {
    Resize(size_half);
  }

  ~Filter1D() {
    delete [] afWeights;
  }

  void Normalize() {
    // Make sure the sum of the filter weights is 1
    RealNum total = 0.0;
    for (Integer ix=-size; ix<=size; ix++)
      total += afWeights[ix+size];
    for (Integer ix=-size; ix<=size; ix++)
      afWeights[ix+size] /= total;
  }

}; // class Filter1D



#endif //#ifndef _FILTER1D_H
