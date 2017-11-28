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
  float window_sigma[3];     //width of Gaussian distribution, sigma (in voxels)
  float window_sigma_nm[3];  //width of Gaussian distribution, sigma (in nm)
  float window_sigma_cutoff; //width of filter window (in units of sigma)
                             //This is an alternative way to specify window_size
  float window_sigma_cutoff_exp; //exp(-0.5*sigma_cutoff^2)
                             
  bool use_thresholds;
  bool use_dual_thresholds;
  float out_threshold_01_a;
  float out_threshold_01_b;
  float out_threshold_10_a;
  float out_threshold_10_b;

  // not used (yet):
  float missing_wedge_min[2];  //range of angles sampled (around x and y axis)
  float missing_wedge_max[2];  //range of angles sampled (around x and y axis)
  Settings();
  void ParseArgs(int argc, char **argv);
 private:
  void ParseArgs(vector<string>& vArgs);
  void ConvertArgvToVectorOfStrings(int argc, 
                                    char **argv, 
                                    vector<string>& dest);
};


#endif //#ifndef _SETTINGS_H
