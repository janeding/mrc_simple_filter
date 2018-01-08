#ifndef _SETTINGS_H
#define _SETTINGS_H



#include<vector>
#include<string>
using namespace std;

// Generally speaking, the filter used is of the form of a
// "difference of gaussians":
//
//   h(x,y,z) = h_a(x,y,z) - h_b(x,y,z)
//
// Both h_a() and h_b() are generalized normal distributions (Gaussians)
// https://en.wikipedia.org/wiki/Generalized_normal_distribution
//
// If b_x>0, b_y>0, b_z>0, then:
//
//   h_a(x,y,z) = A*exp(-r_a^m)
//          r_a = sqrt((x/a_x)*2 + (y/a_y)*2 + (z/a_z)*2)
//   h_b(x,y,z) = B*exp(-r_b^n)
//          r_b = sqrt((x/b_x)*2 + (y/b_y)*2 + (z/b_z)*2)
//
// The "A" and "B" parameters are chosen so that sum_{x,y,z} h(x,y,z) = 0
// (In the special case below where h_b=0, then sum_{x,y,z} h(x,y,z) = 1)
// You may wish to use a "difference of gaussians" approach in some
// directions, and ordinary gaussians" in other directions.
//
// If b_z < 0, then both h_a and h_b become "seperable" filters, IE
//             functions of x,y, multiplied by a Gaussian in the Z direction:
//   h_a(x,y,z) = A * exp(-r_a^m)   *  C * exp(-(1/2)z^2/a_z)
//          r_a = sqrt((x/a_x)*2 + (y/a_y)*2
//   h_b(x,y,z) = B * exp(-r_b^n)   *  C * exp(-(1/2)z^2/a_z)
//          r_b = sqrt((x/b_x)*2 + (y/b_y)*2
// Seperable filters are faster to compute.
// (In this case, "C" is chosen so that C * sum_z exp(-(1/2)z^2/a_z) = 1)
//
// Likewise for each dimension where the corresponding "b" parameter is < 0,
// we use a "separable" filter which is Gaussian in that dimension.
// The second "B" term is only calculated in dimensions for which t>0.
// In this way, the function above is general.
// Thus, a simple 3-D Gaussian (whose principle axes in the x,y,z directions)
// can be implemented setting b_x = b_y = b_z = -1.0.
//   h(x,y,z) =
//     A * exp(-(1/2)(x/a_x)^2) * exp(-(1/2)(y/a_y)^2) * exp(-(1/2)(z/a_z)^2)
// "Seperable" filters like this are faster to compute,
// requiring time  O( Nx * Ny * Nz * max(a_x, a_y, a_z) )
//     instead of  O( Nx * Ny * Nz * a_x * a_y * a_z) )
//                 (where Nx, Ny, Nz denote the source image size in voxels)
//
// Special case:
//       In the specific case where the exponents m = n = 2, this filter
//       can be implemented as the sum of two terms, each 
//       of which are seperable in x,y,z and can be computed quickly:
//   h(x,y,z) = h_a(x,y,z) - h_b(x,y,z)
//   h_a(x,y,z) =
//     A * exp(-(1/2)(x/a_x)^2) * exp(-(1/2)(y/a_y)^2) * exp(-(1/2)(z/a_z)^2)
//   h_b(x,y,z) =
//     B * exp(-(1/2)(x/b_x)^2) * exp(-(1/2)(y/b_y)^2) * exp(-(1/2)(z/b_z)^2)


class Settings {
 public:
  string in_file_name;
  bool in_rescale01;
  bool invert_output;
  string out_file_name;
  string mask_file_name;
  int mask_select;
  bool use_mask_select;
  int mask_out;
  bool use_mask_out;

  float voxel_width; //width of each voxel in nm (assumed to be same for x,y,z)
  bool  voxel_width_divide_by_10 = false;

  float width_a[3];    //"a" parameter in formula above
  float width_b[3];    //"b" parameter in formula above
  float m_exp;         //"m" parameter in formula above
  float n_exp;         //"n" parameter in formula above

  typedef enum eFilterType {
    NONE,    //(the user has not selected a filter)
    GAUSS,   //3D Gaussian filter
    DOG,     //3D difference of Gaussians (fast, separable)
    LOG,     //3D Laplacian-of-Gaussian (fast, separable)
    DOGGEN,  //3D Difference-of-Gaussians filter (general, arbitrary exponents)
    DOGXYGEN //2D Difference-of-Gaussians (in XY direction, arbitrary exponents)
             //and a 1D Gaussian filter (in the Z direction)
  } FilterType; 
  
  FilterType filter_type;

  float window_threshold;
  float window_ratio;
  float window_halfwidth[3];

  bool use_thresholds;
  bool use_dual_thresholds;
  float out_threshold_01_a;
  float out_threshold_01_b;
  float out_threshold_10_a;
  float out_threshold_10_b;

  // not used (yet):
  //float missing_wedge_min[2];  //range of angles sampled (around x and y axis)
  //float missing_wedge_max[2];  //range of angles sampled (around x and y axis)

  Settings();
  void ParseArgs(int argc, char **argv);

 private:
  void ParseArgs(vector<string>& vArgs);
  void ConvertArgvToVectorOfStrings(int argc, 
                                    char **argv, 
                                    vector<string>& dest);
};


#endif //#ifndef _SETTINGS_H
