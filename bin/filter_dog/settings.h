#ifndef _SETTINGS_H
#define _SETTINGS_H



#include<vector>
#include<string>
using namespace std;



class Settings {
 public:
  string in_file_name;
  bool in_rescale01;
  string out_file_name;
  string mask_file_name;
  int mask_select;
  bool use_mask_select;
  int mask_out;
  bool use_mask_out;
  int window_size[3];        //width of filter window (in voxels)

    // Filter used:
  //   h(x,y,z) = h_xy(r) * h_z(z)
  //
  // In the XY plane, the filter used is:
  //   h_xy(r) = A*exp(-(r/s)^m) - B*exp(-(r/t)^n)
  // (Here, "r" refers to disance from the Z axis)
  // (and "A" and "B" are chosen so that h_xy(r=0) = A+B = 1, and
  //  the sum of h_xy(r) for all pixels in the filter range is zero.)
  // Along the Z direction, the filter used is:
  //   h_z(r) = C*exp(-0.5*(z/u)^2)
  // (Where "C" so that the (discrete) sum of h_z(z) values along z is 1)

  float voxel_width; //width of each voxel in nm (assumed to be same for x,y,z)
                             
  float width_s;       //"s" parameter in formula above (in XY plane)
  float width_t;       //"t" parameter in formula above (in XY plane)
  float width_u;       //"u" parameter in formula above (along Zaxis)
  float exponent_m;    //"m" parameter in formula above (in XY plane)
  float exponent_n;    //"n" parameter in formula above (in XY plane)

  float window_sigma_cutoff_xy; //width of filter window (in units of width_t)
  float window_sigma_cutoff_xy_exp; //exp(-(windows_sigma_cutoff_xy)^n)

  float window_sigma_cutoff_z; //width of filter window (in units of width_u)
  float window_sigma_cutoff_z_exp; //exp(-0.5*windows_sigma_cutoff_z^2)

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
