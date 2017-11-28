#include <cmath>
//#include <iostream>
using namespace std;

#include "settings.h"


Settings::Settings() {
  // Default settings
  in_file_name = "";
  in_rescale01 = true;
  out_file_name = "";
  mask_file_name = "";
  mask_select = 1;
  use_mask_select = false;
  mask_out = 0.0;
  use_mask_out = false;
  voxel_width = 0.0;    // How many nm per voxel? (if 0 then read from tomogram)
  window_size[0] = 3;   // A box with (1+2*3)^3 voxels (default, overrideable)
  window_size[1] = 3;   // A box with (1+2*3)^3 voxels
  window_size[2] = 3;   // A box with (1+2*3)^3 voxels

  // Filter used:
  //   h(x,y,z) = h_xy(r) * h_z(z)
  //
  // In the XY plane, the filter used is:
  //
  //   h_xy(r) = A*exp(-(r/s)^m) - B*exp(-(r/t)^n)
  //
  // (Here, "r" refers to disance from the Z axis)
  // (and "A" and "B" are chosen so that h_xy(r=0) = A+B = 1, and
  //  the sum of h_xy(r) for all pixels in the filter range is zero.)
  // Along the Z direction, the filter used is:
  //
  //   h_z(r) = C*exp(-0.5*(z/u)^2)
  //
  // (Where "C" is chosen so the (discrete) sum-over-z of  h_z(z) is 1)

  width_s = 4.0;     //"s" parameter in formula above (in XY plane in nm)
  width_t = 5.0;     //"t" parameter in formula above (in XY plane in nm)
  width_u = 2.5;     //"u" parameter in formula above (along Zaxis in nm)
  exponent_m = 2.0;  //"m" parameter in formula above (in XY plane)
  exponent_n = 2.0;  //"n" parameter in formula above (in XY plane)

  //window_sigma_cutoff_exp_z;
  float window_sigma_cutoff_z=2.0; //By default, when averaging/filtering consider
                             // nearby voxels up to a distance of 2.0*sigma away
                             // Throw away all pixels further than this distance
                             // (even if they lie within the window box).
                             // This setting overrides the window_size[] array.
  use_thresholds = false;
  use_dual_thresholds = false;
  out_threshold_01_a = 1.0;
  out_threshold_01_b = 1.0;
  out_threshold_10_a = 1.0;
  out_threshold_10_b = 1.0;
  //missing_wedge_min[0]=-90.0; // By default, the "missing wedge" includes all
  //missing_wedge_max[0]=90.0;  // orientations around the Y axis which lie
  //missing_wedge_min[1]=-30.0; // between -30 and +30 degrees (relative to the
  //missing_wedge_max[1]=+30.0; // Z axis), independent of X-axis orientation.
}



void Settings::ParseArgs(int argc, char **argv) {
  vector<string> vArgs;
  ConvertArgvToVectorOfStrings(argc, argv, vArgs);
  ParseArgs(vArgs);
}



void Settings::ConvertArgvToVectorOfStrings(int argc,
                                            char **argv,
                                            vector<string>& dest)
{
  dest.resize(argc);
  for (int i=0; i < argc; ++i)
    dest[i].assign(argv[i]);
}



void
Settings::ParseArgs(vector<string>& vArgs)
{
  for (int i=1; i < vArgs.size(); ++i)
  {

    int num_arguments_deleted = 0;

    if (vArgs[i] == "-in")
    {
      if ((i+1 >= vArgs.size()) || (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
        throw "Error: The " + vArgs[i] + 
          " argument must be followed by a file name.\n";
      in_file_name = vArgs[i+1];

      num_arguments_deleted = 2;

    } // if (vArgs[i] == "-in")

    else if (vArgs[i] == "-out")
    {
      if ((i+1 >= vArgs.size()) || (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
        throw "Error: The " + vArgs[i] + 
          " argument must be followed by a file name.\n";
      out_file_name = vArgs[i+1];

      num_arguments_deleted = 2;

    } // if (vArgs[i] == "-out")

    else if (vArgs[i] == "-mask")
    {
      if ((i+1 >= vArgs.size()) || (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
        throw "Error: The " + vArgs[i] + 
          " argument must be followed by a file name.\n";
      mask_file_name = vArgs[i+1];
      num_arguments_deleted = 2;
    } // if (vArgs[i] == "-mask")

    else if (vArgs[i] == "-mask-select")
    {
      if ((i+1 >= vArgs.size()) || (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
        throw "Error: The " + vArgs[i] + 
          " argument must be followed by an integer.\n";
      use_mask_select = true;
      mask_select = stoi(vArgs[i+1]);
      num_arguments_deleted = 2;
    } // if (vArgs[i] == "-mask-select")

    else if (vArgs[i] == "-mask-out")
    {
      if ((i+1 >= vArgs.size()) || (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
        throw "Error: The " + vArgs[i] + 
          " argument must be followed by a number.\n";
      use_mask_out = true;
      mask_out = stof(vArgs[i+1]);
      num_arguments_deleted = 2;
    } // if (vArgs[i] == "-mask-out")

    else if (vArgs[i] == "-w")
    {
      if ((i+1 >= vArgs.size()) || (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
        throw "Error: The " + vArgs[i] + 
          " argument must be followed by voxel width.\n";
      voxel_width = stof(vArgs[i+1]);
      num_arguments_deleted = 2;
    } // if (vArgs[i] == "-mask")

    else if (vArgs[i] == "-norescale")
    {
      in_rescale01 = false;

      num_arguments_deleted = 1;

    } // if (vArgs[i] == "-norescale")

    else if ((vArgs[i] == "-exponent") || (vArgs[i] == "-exponents"))
    {
      if ((i+2 >= vArgs.size()) ||
          (vArgs[i+1] == "") || (vArgs[i+1][0] == '-') ||
          (vArgs[i+2] == "") || (vArgs[i+2][0] == '-'))
        throw "Error: The " + vArgs[i] + 
          " argument must be followed by two positive numbers.\n";
      exponent_m = stof(vArgs[i+1]);
      exponent_n = stof(vArgs[i+2]);
      num_arguments_deleted = 3;
    } //if (vArgs[i] == "-width-z")

    else if (vArgs[i] == "-width-xy")
    {
      if ((i+2 >= vArgs.size()) ||
          (vArgs[i+1] == "") || (vArgs[i+1][0] == '-') ||
          (vArgs[i+2] == "") || (vArgs[i+2][0] == '-'))
        throw string("Error: The " + vArgs[i] + 
                     " argument must be followed by 2 positive numbers.\n");
      width_s = stof(vArgs[i+1]);
      width_t = stof(vArgs[i+2]);
      if (width_t <= width_s)
        throw string("Error: The two arguments following " + vArgs[i] + " must be\n"
                     " increasing.  (Ie., the 2nd argument must be > 1st argument.)\n");
      num_arguments_deleted = 3;
    } //if (vArgs[i] == "-width-xy")

    else if (vArgs[i] == "-width-z")
    {
      if ((i+1 >= vArgs.size()) ||
          (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
        throw "Error: The " + vArgs[i] + 
          " argument must be followed by a positive number.\n";
      width_u = stof(vArgs[i+1]);
      num_arguments_deleted = 2;
    } //if (vArgs[i] == "-width-z")

    else if (vArgs[i] == "-cutoff-xy")
    {
      if ((i+1 >= vArgs.size()) || (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
        throw "Error: The " + vArgs[i] + 
          " argument must be followed by a number.\n";
      window_sigma_cutoff_xy = stof(vArgs[i+1]);
      window_sigma_cutoff_xy_exp = exp(-pow(window_sigma_cutoff_xy,
                                            exponent_n));
      num_arguments_deleted = 2;
    } // if (vArgs[i] == "-cutoff-xy")

    else if (vArgs[i] == "-cutoff-z")
    {
      if ((i+1 >= vArgs.size()) || (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
        throw "Error: The " + vArgs[i] + 
          " argument must be followed by a number.\n";
      window_sigma_cutoff_z = stof(vArgs[i+1]);
      window_sigma_cutoff_z_exp = exp(-0.5 * window_sigma_cutoff_z *
                                             window_sigma_cutoff_z);
      num_arguments_deleted = 2;
    } // if (vArgs[i] == "-cutoff-z")

    else if (vArgs[i] == "-window") {
      if ((i+1 >= vArgs.size()) || (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
        throw "Error: The " + vArgs[i] + 
          " argument must be followed by either 1 or 3 positive integers.\n";
      window_size[0] = stoi(vArgs[i+1]);
      num_arguments_deleted = 2;
      if ((i+4 >= vArgs.size()) && 
          (vArgs[i+2][0] != '-') && (vArgs[i+2] != "") &&
          (vArgs[i+3][0] != '-') && (vArgs[i+3] != "")) {
        window_size[1] = stoi(vArgs[i+2]);
        window_size[2] = stoi(vArgs[i+3]);
        num_arguments_deleted = 4;
      }
      else {
        window_size[1] = window_size[0];
        window_size[2] = window_size[0];
      }
    } //if (vArgs[i] == "-window")

    else if (vArgs[i] == "-thresh") {
      if (i+2 >= vArgs.size())
        throw "Error: The " + vArgs[i] + 
          " argument must be followed by 1 number.\n";
      use_thresholds = true;
      use_dual_thresholds = false;
      out_threshold_01_a = stof(vArgs[i+1]);
      out_threshold_01_b = out_threshold_01_a;
    }
    else if (vArgs[i] == "-thresh2") {
      if (i+3 >= vArgs.size())
        throw "Error: The " + vArgs[i] + 
          " argument must be followed by 2 numbers.\n";
      use_thresholds = true;
      use_dual_thresholds = false;
      out_threshold_01_a = stof(vArgs[i+1]);
      out_threshold_01_b = stof(vArgs[i+2]);
    }
    else if (vArgs[i] == "-thresh4") {
      if (i+5 >= vArgs.size())
        throw "Error: The " + vArgs[i] + 
          " argument must be followed by 4 numbers.\n";
      use_thresholds = true;
      use_dual_thresholds = true;
      out_threshold_01_a = stof(vArgs[i+1]);
      out_threshold_01_b = stof(vArgs[i+2]);
      out_threshold_10_a = stof(vArgs[i+3]);
      out_threshold_10_b = stof(vArgs[i+4]);
    }


    // NOT USED YET:
    //else if (vArgs[i] == "-xwedge") {
    //  if ((i+1 >= vArgs.size()) || (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
    //    throw "Error: The " + vArgs[i] + 
    //      " argument must be followed by 1 or 3 positive numbers.\n";
    //  missing_wedge_min[0] = stof(vArgs[i+1]);
    //  num_arguments_deleted = 2;
    //  if ((i+4 >= vArgs.size()) && 
    //      (vArgs[i+2][0] != '-') && (vArgs[i+2] != "")) {
    //    missing_wedge_max[0] = stof(vArgs[i+2]);
    //    num_arguments_deleted = 3;
    //  }
    //  else {
    //    missing_wedge_max[0] = missing_wedge_min[0];
    //    missing_wedge_min[0] = -missing_wedge_min[0];
    //  }
    //} //if (vArgs[i] == "-xwedge")
    //
    //else if ((vArgs[i] == "-ywedge") || (vArgs[i] == "-wedge")) {
    //  if ((i+1 >= vArgs.size()) || (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
    //    throw "Error: The " + vArgs[i] + 
    //      " argument must be followed by 1 or 3 positive numbers.\n";
    //  missing_wedge_min[1] = stof(vArgs[i+1]);
    //  num_arguments_deleted = 2;
    //  if ((i+4 >= vArgs.size()) && 
    //      (vArgs[i+2][0] != '-') && (vArgs[i+2] != "")) {
    //    missing_wedge_max[1] = stof(vArgs[i+2]);
    //    num_arguments_deleted = 3;
    //  }
    //  else {
    //    missing_wedge_max[1] = missing_wedge_min[1];
    //    missing_wedge_min[1] = -missing_wedge_min[1];
    //  }
    //} //if (vArgs[i] == "-ywedge")



    //Delete all the arguments we have processed in this iteration (if any)
    if (num_arguments_deleted > 0)
    {
      vArgs.erase(vArgs.begin()+i, vArgs.begin()+i+num_arguments_deleted);
      --i; //rewind by one so that when we invoke ++i
      //at the end of this loop, i's value will not change.
      //This will point us to the next un-read argument.
    }
    
  } // loop over arguments "for (int i=1; i < vArgs.size(); ++i)"

} // Settings::ParseArgs()



