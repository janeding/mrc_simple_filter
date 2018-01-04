#include <cmath>
//#include <iostream>
using namespace std;

#include "settings.h"


template<class RealNum>
inline RealNum MAX(RealNum x, RealNum y) { return ((x<y) ? y : x); }


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
  voxel_width_divide_by_10 = false;
  filter_type = GAUSS;

  // The code is a little bit confusing because I tried to cram as much
  // functionality into the smallest number of parameters possible:
  // The "difference of gaussians" (dog) filter used here requires
  //  8 parameters: a_x, a_y, a_z, b_x, b_y, b_z, m, n
  // For details, see the comment at the beginning of "settings.h".
  //
  // Default values of a_x, a_y, a_z, b_x, b_y, b_z, m, n are:

  width_a[0] = 12.0; // a_x = 12.0  gaussian width in x direction
  width_a[1] = 12.0; // a_y = 12.0  gaussian width in y direction
  width_a[2] = 24.0; // a_z = 24.0  gaussian width in z direction
  width_b[0] = 15.0; // b_x = 15.0  gaussian width in x direction
  width_b[1] = 15.0; // b_y = 15.0  gaussian width in y direction
  width_b[2] = -1.0; // b_z = -1.0  gaussian width in z direction (<0 disables)
  m_exp = 2.0;       // exponent in generalized Gaussian formula (width a)
  n_exp = 2.0;       // exponent in generalized Gaussian formula (width b)

  filter_halfwidth[0] = -1; // width of the filter used (in voxels)
  filter_halfwidth[1] = -1; // (This is the size of the domain of the function
  filter_halfwidth[2] = -1; //  which will be convolved with the image.
                              //  "-1" means unspecified.)

  window_cutoff_ratio = 3.0;
                             //By default, when averaging/filtering consider
                             // nearby voxels up to a distance of 3.0*sigma away
                             // Throw away all pixels further than this distance
                             // (even if they lie within the window box).
                             //This only occurs if the filter_halfwidth was
                             //not otherwise specified manually by the user.

  //float window_sigma_dog = 2.0;   //In directions where a diff-of-gauss is used
  //                                //use this as the window size (in units of
  //                                //the "b" parameter in that direction).
  //float window_sigma_gauss = 2.0; //In directions where a Gaussian is used,
  //                                //use this as the window size (in units of
  //                                //the "a" (sigma) parameter in that direction)
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
  bool exponents_set_by_user = false;
  float log_width[3];
  float log_delta_t_over_t = 0.01;
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


    else if (vArgs[i] == "-ang-to-nm")
    {
      voxel_width_divide_by_10 = true;
      num_arguments_deleted = 1;
    } // if (vArgs[i] == "-mask")


    else if (vArgs[i] == "-norescale")
    {
      in_rescale01 = false;

      num_arguments_deleted = 1;

    } // if (vArgs[i] == "-norescale")


    else if (vArgs[i] == "-gauss-aniso")
    {
      if ((i+3 >= vArgs.size()) ||
          (vArgs[i+1] == "") || (vArgs[i+1][0] == '-') ||
          (vArgs[i+2] == "") || (vArgs[i+2][0] == '-') ||
          (vArgs[i+3] == "") || (vArgs[i+3][0] == '-'))
        throw string("Error: The " + vArgs[i] + 
                     " argument must be followed by 3 positive numbers:\n"
                     " a_xy  b_xy  a_z\n"
                     " (I.E., the \"A\" and \"-B\" Gaussian widths in the XY plane,\n"
                     "  followed by the Gaussian width in the Z direction.)\n");
      width_a[0] = stof(vArgs[i+1]);
      width_a[1] = stof(vArgs[i+2]);
      width_a[2] = stof(vArgs[i+3]);
      width_b[0] = -1.0;
      width_b[1] = -1.0;
      width_b[2] = -1.0;
      filter_type = GAUSS;
      num_arguments_deleted = 4;
    } //if (vArgs[i] == "-gauss-aniso")

    else if ((vArgs[i] == "-gauss") || (vArgs[i] == "-gauss-iso"))
    {
      if ((i+1 >= vArgs.size()) ||
          (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
        throw string("Error: The " + vArgs[i] + 
                     " argument must be followed by 3 positive numbers:\n"
                     " a_xy  b_xy  a_z\n"
                     " (I.E., the \"A\" and \"-B\" Gaussian widths in the XY plane,\n"
                     "  followed by the Gaussian width in the Z direction.)\n");
      width_a[0] = stof(vArgs[i+1]);
      width_a[1] = width_a[0];
      width_a[2] = width_a[0];
      width_b[0] = -1.0;
      width_b[1] = -1.0;
      width_b[2] = -1.0;
      filter_type = GAUSS;
      num_arguments_deleted = 2;
    } //if (vArgs[i] == "-gauss")


    else if ((vArgs[i] == "-dog-aniso") || (vArgs[i] == "-dogxyz-aniso"))
    {
      if ((i+6 >= vArgs.size()) ||
          (vArgs[i+1] == "") || (vArgs[i+1][0] == '-') ||
          (vArgs[i+2] == "") || (vArgs[i+2][0] == '-') ||
          (vArgs[i+3] == "") || (vArgs[i+3][0] == '-') ||
          (vArgs[i+4] == "") || (vArgs[i+4][0] == '-') ||
          (vArgs[i+5] == "") || (vArgs[i+5][0] == '-') ||
          (vArgs[i+6] == "") || (vArgs[i+6][0] == '-'))
        throw string("Error: The " + vArgs[i] + 
                     " argument must be followed by 3 positive numbers.\n");
      width_a[0] = stof(vArgs[i+1]);
      width_a[1] = stof(vArgs[i+2]);
      width_a[2] = stof(vArgs[i+3]);
      width_b[0] = stof(vArgs[i+4]);
      width_b[1] = stof(vArgs[i+5]);
      width_b[2] = stof(vArgs[i+6]);
      //if (width_b[0] <= width_a[0])
      //  throw string("Error: The two arguments following " + vArgs[i] + " must be\n"
      //               " increasing.  (Ie., the 2nd argument must be > 1st argument.)\n");
      filter_type = DOGGEN;
      num_arguments_deleted = 7;
    } //if (vArgs[i] == "-dog")


    else if ((vArgs[i] == "-dog") || (vArgs[i] == "-dog-iso"))
    {
      if ((i+2 >= vArgs.size()) ||
          (vArgs[i+1] == "") || (vArgs[i+1][0] == '-') ||
          (vArgs[i+2] == "") || (vArgs[i+2][0] == '-'))
        throw string("Error: The " + vArgs[i] + 
                     " argument must be followed by 3 positive numbers.\n");
      width_a[0] = stof(vArgs[i+1]);
      width_a[1] = width_a[0];
      width_a[2] = width_a[0];
      width_b[0] = stof(vArgs[i+2]);
      width_b[1] = width_b[0];
      width_b[2] = width_b[0];
      //if (width_b[0] <= width_a[0])
      //  throw string("Error: The two arguments following " + vArgs[i] + " must be\n"
      //               " increasing.  (Ie., the 2nd argument must be > 1st argument.)\n");
      filter_type = DOGGEN;
      num_arguments_deleted = 3;
    } //if (vArgs[i] == "-dog")

    else if (vArgs[i] == "-dogxy-aniso")
    {
      if ((i+5 >= vArgs.size()) ||
          (vArgs[i+1] == "") || (vArgs[i+1][0] == '-') ||
          (vArgs[i+2] == "") || (vArgs[i+2][0] == '-') ||
          (vArgs[i+3] == "") || (vArgs[i+3][0] == '-') ||
          (vArgs[i+4] == "") || (vArgs[i+4][0] == '-') ||
          (vArgs[i+5] == "") || (vArgs[i+5][0] == '-'))
        throw string("Error: The " + vArgs[i] + 
                     " argument must be followed by 5 positive numbers:\n"
                     " a_x  a_y  b_x  b_y  a_z\n"
                     " (I.E., the \"A\" and \"-B\" Gaussian widths in X and Y directions,\n"
                     "  followed by the Gaussian width in the Z direction.)\n");
      width_a[0] = stof(vArgs[i+1]);
      width_a[1] = stof(vArgs[i+2]);
      width_b[0] = stof(vArgs[i+3]);
      width_b[1] = stof(vArgs[i+4]);
      //The "-dogxy" filter is a difference of Gaussians in the X,Y directions
      //multiplied by a Gaussian in the Z direction.
      //By default, this Gaussian has width 0.
      //(In other words, by default it has no blurring in the Z direction.)
      //Users can override this using the "-gauss-z" argument.
      width_a[2] = stof(vArgs[i+5]);
      width_b[2] = -1.0;   //(This disables the second (negative) Gaussian)
      //if (width_b[0] <= width_a[0])
      //  throw string("Error: The two arguments following " + vArgs[i] + " must be\n"
      //               " increasing.  (Ie., the 2nd argument must be > 1st argument.)\n");
      filter_type = DOGXYGEN;
      num_arguments_deleted = 6;
    } //if (vArgs[i] == "-dogxy-aniso")


    else if ((vArgs[i] == "-dogxy") || (vArgs[i] == "-dogxy-iso"))
    {
      if ((i+5 >= vArgs.size()) ||
          (vArgs[i+1] == "") || (vArgs[i+1][0] == '-') ||
          (vArgs[i+2] == "") || (vArgs[i+2][0] == '-') ||
          (vArgs[i+3] == "") || (vArgs[i+3][0] == '-'))
        throw string("Error: The " + vArgs[i] + 
                     " argument must be followed by 3 positive numbers:\n"
                     " a_xy  b_xy  a_z\n"
                     " (I.E., the \"A\" and \"-B\" Gaussian widths in the XY plane,\n"
                     "  followed by the Gaussian width in the Z direction.)\n");
      width_a[0] = stof(vArgs[i+1]);
      width_a[1] = width_a[0];
      width_b[0] = stof(vArgs[i+2]);
      width_b[1] = width_b[0];
      //The "-dogxy" filter is a difference of Gaussians in the X,Y directions
      //multiplied by a Gaussian in the Z direction.
      //By default, this Gaussian has width 0.
      //(In other words, by default it has no blurring in the Z direction.)
      //Users can override this using the "-gauss-z" argument.
      width_a[2] = stof(vArgs[i+3]);
      width_b[2] = -1.0;   //(This disables the second (negative) Gaussian)
      //if (width_b[0] <= width_a[0])
      //  throw string("Error: The two arguments following " + vArgs[i] + " must be\n"
      //               " increasing.  (Ie., the 2nd argument must be > 1st argument.)\n");
      filter_type = DOGXYGEN;
      num_arguments_deleted = 4;
    } //if (vArgs[i] == "-dogxy")


    else if ((vArgs[i] == "-log-aniso") || (vArgs[i] == "-blob-aniso"))
    {
      if ((i+3 >= vArgs.size()) ||
          (vArgs[i+1] == "") || (vArgs[i+1][0] == '-') ||
          (vArgs[i+2] == "") || (vArgs[i+2][0] == '-') ||
          (vArgs[i+3] == "") || (vArgs[i+3][0] == '-'))
        throw string("Error: The " + vArgs[i] + 
                     " argument must be followed by 3 positive numbers.\n");
      log_width[0] = stof(vArgs[i+1]);
      log_width[1] = stof(vArgs[i+2]);
      log_width[2] = stof(vArgs[i+3]);
      filter_type = LOG;
      num_arguments_deleted = 4;
    } //if ((vArgs[i] == "-log-aniso") || (vArgs[i] == "-blob-aniso"))


    else if ((vArgs[i] == "-log") ||
             (vArgs[i] == "-log-iso") ||
             (vArgs[i] == "-blob") ||
             (vArgs[i] == "-blob-iso"))
    {
      if ((i+1 >= vArgs.size()) ||
          (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
        throw string("Error: The " + vArgs[i] + 
                     " argument must be followed by 3 positive numbers.\n");
      log_width[0] = stof(vArgs[i+1]);
      log_width[1] = log_width[0];
      log_width[2] = log_width[0];
      filter_type = LOG;
      num_arguments_deleted = 2;
    } //if ((vArgs[i] == "-blob-aniso") || (vArgs[i] == "-log-aniso"))


    else if ((vArgs[i] == "-blob-delta") || (vArgs[i] == "-log-delta")) {
      if ((i+1 >= vArgs.size()) ||
          (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
        throw string("Error: The " + vArgs[i] + 
                     " argument must be followed by 1 positive number.\n");
      log_delta_t_over_t = stof(vArgs[i+1]);
    }


    else if (vArgs[i] == "-cutoff")
    {
      if ((i+3 >= vArgs.size()) || (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
        throw "Error: The " + vArgs[i] + 
          " argument must be followed by a number.\n";
      window_cutoff_ratio = stof(vArgs[i+1]);
      //window_cutoff_ratio_exp = exp(-pow(window_cutoff_ratio,
      //                                   n_exp));
      num_arguments_deleted = 2;
    } // if (vArgs[i] == "-cutoff")


    //else if (vArgs[i] == "-cutoff-dog")
    //{
    //  if ((i+1 >= vArgs.size()) || (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
    //    throw "Error: The " + vArgs[i] + 
    //      " argument must be followed by a number.\n";
    //  window_cutoff_ratio_dog = stof(vArgs[i+1]);
    //  //window_cutoff_ratio_dog = exp(-pow(window_sigmma_dog, n_exp));
    //  num_arguments_deleted = 2;
    //} // if (vArgs[i] == "-cutoff-dog")


    //else if (vArgs[i] == "-cutoff-gauss")
    //{
    //  if ((i+1 >= vArgs.size()) || (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
    //    throw "Error: The " + vArgs[i] + 
    //      " argument must be followed by a number.\n";
    //  window_cutoff_ratio_gauss = stof(vArgs[i+1]);
    //  //window_cutoff_ratio_gauss_exp = exp(-0.5 * window_sigma_gauss *
    //  //                                           window_sigma_gauss);
    //  num_arguments_deleted = 2;
    //} // if (vArgs[i] == "-cutoff-z")


    else if ((vArgs[i] == "-dog-exponents") ||
             (vArgs[i] == "-exponents"))
    {
      if ((i+2 >= vArgs.size()) ||
          (vArgs[i+1] == "") || (vArgs[i+1][0] == '-') ||
          (vArgs[i+2] == "") || (vArgs[i+2][0] == '-'))
        throw "Error: The " + vArgs[i] + 
          " argument must be followed by two positive numbers.\n";
      m_exp = stof(vArgs[i+1]);
      n_exp = stof(vArgs[i+2]);
      exponents_set_by_user = true;
      num_arguments_deleted = 3;
    } //if (vArgs[i] == "-dog-exponents")


    else if ((vArgs[i] == "-gauss-exponent") ||
             (vArgs[i] == "-exponent"))
    {
      if ((i+1 >= vArgs.size()) ||
          (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
        throw "Error: The " + vArgs[i] + 
          " argument must be followed by two positive numbers.\n";
      m_exp = stof(vArgs[i+1]);
      n_exp = m_exp;
      exponents_set_by_user = true;
      num_arguments_deleted = 2;
    } //if (vArgs[i] == "-gauss-exponent")


    else if (vArgs[i] == "-window") {
      if ((i+1 >= vArgs.size()) || (vArgs[i+1] == "") || (vArgs[i+1][0] == '-'))
        throw "Error: The " + vArgs[i] + 
          " argument must be followed by either 1 or 3 positive integers.\n";
      filter_halfwidth[0] = stoi(vArgs[i+1]);
      num_arguments_deleted = 2;
      if ((i+4 >= vArgs.size()) && 
          (vArgs[i+2][0] != '-') && (vArgs[i+2] != "") &&
          (vArgs[i+3][0] != '-') && (vArgs[i+3] != "")) {
        filter_halfwidth[1] = stoi(vArgs[i+2]);
        filter_halfwidth[2] = stoi(vArgs[i+3]);
        num_arguments_deleted = 4;
      }
      else {
        filter_halfwidth[1] = filter_halfwidth[0];
        filter_halfwidth[2] = filter_halfwidth[0];
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



  // ----------
  if (in_file_name.size() == 0) {
    throw "Error: You must specify the name of the tomogram you want to read\n"
          "       using the \"-in SOURCE_FILE.mrc\" argument\n"
          "       (\".rec\" files are also supported.)\n";
  }
  if (out_file_name.size() == 0) {
    throw "Error: You must specify the name of the tomogram you want to read\n"
          "       using the \"-out DESTINATION_FILE.mrc\" argument\n"
          "       (\".rec\" files are also supported.)\n";
  }

  if (filter_type == LOG) {
    // "-log" approximates to the "Laplacian of a Gaussian" ("LOG") filter
    // with the difference of two Gaussians ("DOG") filter.
    // (The two Gaussians have widths which are slightly above and
    //  slightly below the width parameters specified by the user.)
    // For background details see:
    // https://en.wikipedia.org/wiki/Blob_detection
    // https://en.wikipedia.org/wiki/Difference_of_Gaussians
    // https://en.wikipedia.org/wiki/Mexican_hat_wavelet
    width_a[0] = log_width[0] * (1.0 - 0.5*log_delta_t_over_t);
    width_a[1] = log_width[1] * (1.0 - 0.5*log_delta_t_over_t);
    width_a[2] = log_width[2] * (1.0 - 0.5*log_delta_t_over_t);
    width_b[0] = log_width[0] * (1.0 + 0.5*log_delta_t_over_t);
    width_b[1] = log_width[1] * (1.0 + 0.5*log_delta_t_over_t);
    width_b[2] = log_width[2] * (1.0 + 0.5*log_delta_t_over_t);
    m_exp = 2.0;
    n_exp = 2.0;
    filter_type = DOG;
  }

  // ----------
  if ((filter_type == DOGGEN) &&
      (m_exp == 2.0) &&
      (n_exp == 2.0)) {

    filter_type = DOG;

    if (exponents_set_by_user) {
      //The default settings include a factor of 1/sqrt(2) in the gaussian width
      //The formula for a normal Gaussian distribution (used by default)includes
      // a factor of 1/2.  The formula for the generalized Gaussian does not.
      // If the user manually specified the exponents, then I assume they
      // don't want us to include this factor of 1/sqrt(2), even if they choose
      // the exponents to be 2 (as they would be under default conditions 
      // using normal Gaussians.) Otherwise, they may be surprised why something
      // special happens when the exponents = 2.0 exactly.  Later we will check
      // if the expoents == 2.0, and use a faster algorithm (DOG) which 
      // includes this factor of 1/sqrt(2) automatically.
      // So we compensate for that now:
      for (int d=0; d<3; d++) {
        if (width_a[d] > 0)
          width_a[d] *= sqrt(2.0);
        if (width_b[d] > 0)
          width_b[d] *= sqrt(2.0);
      }
    }
  } //if ((filter_type == DOGGEN) && (m_exp == 2.0) && (m_exp == 2.0))


  // ----------

  // If the user did not specify the width of the filters explicitly,
  // then, determine the filter window size (filter width) from the
  // "window_cutoff_ratio" parameters which are distances expressed
  // as multiples of the widths of the gaussians (a and b parameters)
  for (int d=0; d < 3; d++) {
    if (filter_halfwidth[d] < 0.0)
      filter_halfwidth[d] = (MAX(width_a[d], width_b[d])
                             *
                             window_cutoff_ratio);
      //(NOTE: These units are still in nm (physical units), not voxels.
      // We still need to convert them into voxels. When we do we will use the 
      // "ceil()" function to make sure all of the windows have integer size.)
  }

} // Settings::ParseArgs()



