#ifndef _FILTER1D_H
#define _FILTER1D_H

using namespace std;




template<class RealNum, class Integer>

class Filter1D {
public:
  RealNum *afWeights;
  Integer halfwidth; //distance from the filter center to its edge, in pixels

  // Apply the filter to tomographic data in the "afSource" array
  // Save the results in the "afDest" array.  (A "mask" is optional.)
  // All arrays are 1D and assumed to be the same size
  void Apply(Integer const size_source,
             RealNum *afSource,
             RealNum *afDest,
             RealNum const *afMask = NULL) const
             //bool precompute_mask_times_source = true) const
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

    // The mask should be 1 everywhere we want to consider, and 0 elsewhere.
    // Multiplying the density in the tomogram by the mask removes some of 
    // the voxels from consideration later on when we do the filtering.
    // (Later, we will adjust the weight of the average we compute when we
    //  apply the filter in order to account for the voxels we deleted now.)
    // Precomupting the mask product is faster but modifies the source image.

    //if (afMask && precompute_mask_times_source)
    if (afMask) 
      for (int ix=0; ix<size_source; ix++)
        afSource[ix] *= afMask[ix];

    for (Integer ix=0; ix<size_source; ix++) {

      if ((afMask) && (afMask[ix] == 0.0))
        continue;
          
      RealNum g = 0.0;
      RealNum denominator = 0.0;

      for (Integer jx=-halfwidth; jx<=halfwidth; jx++) {

        if ((ix-jx < 0) || (size_source <= ix-jx))
          continue;

        RealNum delta_g = afWeights[jx+halfwidth] * afSource[ix-jx];

        //if (! precompute_mask_times_source)
        //  delta_g *= afMask[ix-jx];

        g += delta_g;
          // Note: We previously applied the mask by multiplying
          //          aDensity[ix]
          //         by afMask[ix]   (if present)
        if (afMask)
          denominator += afMask[ix-jx] * afWeights[jx+halfwidth];
        else
          denominator += afWeights[jx+halfwidth];
                                          
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

  void Resize(Integer set_halfwidth) {
    if (afWeights)
      delete [] afWeights;
    halfwidth = set_halfwidth;
    Integer array_size = 1 + 2*halfwidth;
    afWeights = new RealNum [array_size];
  }

  Filter1D() {
    halfwidth = -1;
    afWeights = NULL;
  }

  Filter1D(Integer halfwidth) {
    Resize(halfwidth);
  }

  ~Filter1D() {
    delete [] afWeights;
  }

  void Normalize() {
    // Make sure the sum of the filter weights is 1
    RealNum total = 0.0;
    for (Integer ix=-halfwidth; ix<=halfwidth; ix++)
      total += afWeights[ix+halfwidth];
    for (Integer ix=-halfwidth; ix<=halfwidth; ix++)
      afWeights[ix+halfwidth] /= total;
  }

}; // class Filter1D


template<class RealNum, class Integer>
inline Filter1D<RealNum, Integer>&
operator = (const Filter1D<RealNum, Integer>& source) {
  Dealloc(); // (just in case)
  Resize(source.halfwidth); // allocates and initializes af and aaafWeights
  //for(Int ix=-halfwidth; ix<=halfwidth; ix++)
  //  afWeights[ix] = source.afWeights[ix];
  // Use memcpy() instead:
  memcpy(afWeights,
         source.afWeights,
         (1+2*halfwidth[0]) * sizeof(RealNum));
}

#endif //#ifndef _FILTER1D_H
