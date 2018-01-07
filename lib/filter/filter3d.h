#ifndef _FILTER3D_H
#define _FILTER3D_H

#include <ostream>
using namespace std;





template<class RealNum, class Integer>

class Filter3D {
public:
  RealNum *af;
  RealNum ***aaafWeights;
  Integer halfwidth[3]; //num pixels from filter center to edge in x,y,z directions



  // Apply the filter to tomographic data in the "aaafSource" array
  // Save the results in the "aaafDest" array.  (A "mask" is optional.)
  // All arrays are 3D and assumed to be the same size
  void Apply(Integer const size_source[3],
             RealNum ***aaafSource,
             RealNum ***aaafDest,
             RealNum ***aaafMask = NULL,
             //bool precompute_mask_times_source = true,
             ostream *preport_progress = NULL) const
  {

    // Apply the filter to the original tomogram data. (Store in aaafSource)
    //        ___
    //        \
    // g(i) = /__  h(j) * f(i-j)
    //         j
    //
    // where: f(i) is the original density of the tomogram at position ix,iy,iz
    //        h(j) is the filter (smoothing function)
    //       
    //  Note on mask functions:
    //     When summing over j, we ignore contributions from voxels where
    //     mask(i-j) is zero.  We don't count them in the average.
    //     Because h(j) is not necessarily normalized, g(i) is later divided by
    //     the area under the curve h(j)*mask(i-j)    (as a function of j)
    //     

    if (preport_progress)
      *preport_progress << "  progress: processing plane#" << endl;

    // The mask should be 1 everywhere we want to consider, and 0 elsewhere.
    // Multiplying the density in the tomogram by the mask removes some of 
    // the voxels from consideration later on when we do the filtering.
    // (Later, we will adjust the weight of the average we compute when we
    //  apply the filter in order to account for the voxels we deleted now.)
    // Precomupting the mask product is faster but modifies the source image.

    //if (aaafMask && precompute_mask_times_source)
    if (aaafMask) 
      for (int iz=0; iz<size_source[2]; iz++)
        for (int iy=0; iy<size_source[1]; iy++)
          for (int ix=0; ix<size_source[0]; ix++)
            aaafSource[iz][iy][ix] *= aaafMask[iz][iy][ix];

    for (Integer iz=0; iz<size_source[2]; iz++) {

      if (preport_progress)
        *preport_progress << "  " << iz+1 << " / " << size_source[2] << "\n";

      for (Integer iy=0; iy<size_source[1]; iy++) {

        for (Integer ix=0; ix<size_source[0]; ix++) {

          if ((aaafMask) && (aaafMask[iz][iy][ix] == 0.0))
          continue;
          
          RealNum g = 0.0;
          RealNum denominator = 0.0;

          for (Integer jz=-halfwidth[2]; jz<=halfwidth[2]; jz++) {

            if ((iz-jz < 0) || (size_source[2] <= iz-jz))
              continue;

            for (Integer jy=-halfwidth[1]; jy<=halfwidth[1]; jy++) {

              if ((iy-jy < 0) || (size_source[1] <= iy-jy))
                continue;

              for (Integer jx=-halfwidth[0]; jx<=halfwidth[0]; jx++) {

                if ((ix-jx < 0) || (size_source[0] <= ix-jx))
                  continue;

                RealNum delta_g = 
                  aaafWeights[jz+halfwidth[2]]
                             [jy+halfwidth[1]]
                             [jx+halfwidth[0]] *
                   aaafSource[iz-jz][iy-jy][ix-jx];

                //if (! precompute_mask_times_source)
                //  delta_g *= aaafMask[iz-jz][iy-jy][ix-jx];

                g += delta_g;
                       // Note: We previously applied the mask by multiplying
                       //          aaafDensity[iz][iy][ix]
                       //          by aaafMask[iz][iy][ix]   (if present)
                if (aaafMask)
                  denominator += aaafMask[iz-jz][iy-jy][ix-jx] *
                                      aaafWeights[jz+halfwidth[2]]
                                                 [jy+halfwidth[1]]
                                                 [jx+halfwidth[0]];
                else
                  denominator += aaafWeights[jz+halfwidth[2]]
                                            [jy+halfwidth[1]]
                                            [jx+halfwidth[0]];
                // Note: If there were no mask, and if the filter is normalized
                // then denominator=1 always, and we could skip the line above.
              }
            }
          }


          if (denominator > 0.0)
            g /= denominator;
          else
            //Otherwise, this position lies outside the mask region.
            g = 0.0;

          aaafDest[iz][iy][ix] = g;

        }
      }
    }
  } //Filter3D::Apply()


  void Resize(Integer const set_halfwidth[3]) {
    Integer array_size[3];
    if (af && aaafWeights) {
      for(Integer d=0; d < 3; d++)
        array_size[d] = 1 + 2*halfwidth[d];
      Dealloc3D(array_size, &af, &aaafWeights);
    }
    for(Integer d=0; d < 3; d++) {
      halfwidth[d] = set_halfwidth[d];
      array_size[d] = 1 + 2*halfwidth[d];
    }
    Alloc3D(array_size, &af, &aaafWeights);
  }

  Filter3D(Integer const set_halfwidth[3]) {
    af = NULL;
    aaafWeights = NULL;
    Resize(set_halfwidth);
  }

  Filter3D() {
    halfwidth[0] = -1;
    halfwidth[1] = -1;
    halfwidth[2] = -1;
    af = NULL;
    aaafWeights = NULL;
  }

  ~Filter3D() {
    Integer array_size[3];
    for(Integer d=0; d < 3; d++)
      array_size[d] = 1 + 2*halfwidth[d];
    Dealloc3D(array_size, &af, &aaafWeights);
  }

  void Normalize() {
    // Make sure the sum of the filter weights is 1
    RealNum total = 0.0;
    for (Integer iz=-halfwidth[2]; iz<=halfwidth[2]; iz++)
      for (Integer iy=-halfwidth[1]; iy<=halfwidth[1]; iy++)
        for (Integer ix=-halfwidth[0]; ix<=halfwidth[0]; ix++)
          total +=
            aaafWeights[iz+halfwidth[2]][iy+halfwidth[1]][ix+halfwidth[0]];
    for (Integer iz=-halfwidth[2]; iz<=halfwidth[2]; iz++)
      for (Integer iy=-halfwidth[1]; iy<=halfwidth[1]; iy++)
        for (Integer ix=-halfwidth[0]; ix<=halfwidth[0]; ix++)
          aaafWeights[iz+halfwidth[2]][iy+halfwidth[1]][ix+halfwidth[0]]
            /= total;
  }
}; // class Filter3D


template<class RealNum, class Integer>
inline Filter3D<RealNum, Integer>&
operator = (const Filter3D<RealNum, Integer>& source) {
  Dealloc(); // (just in case)
  Resize(source.halfwidth); // allocates and initializes af and aaafWeights
  //for(Int iz=-halfwidth[2]; iz<=halfwidth[2]; iz++)
  //  for(Int iy=-halfwidth[1]; iy<=halfwidth[1]; iy++)
  //    for(Int ix=-halfwidth[0]; ix<=halfwidth[0]; ix++)
  //      aaafWeights[iz][iy][ix] = source.aaafWeights[iz][iy][ix];
  // Use memcpy() instead:
  memcpy(af,
         source.af,
         ((1+2*halfwidth[0]) * (1+2*halfwidth[1]) * (1+2*halfwidth[2]))
         *sizeof(RealNum));
}

#endif //#ifndef _FILTER3D_H
