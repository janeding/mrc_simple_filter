#include <string>
#include <iostream>
using namespace std;
#include "mrc_simple.h"

// (Note: For gcc version 4.8.3, you must compile using: g++ -std=c++11)

string g_program_name = "crop_mrc";

int main(int argc, char **argv) {
  if (argc <= 8) {
    cerr << "This program expects 8 arguments.  Syntax:\n"
	 << "  " << g_program_name << " input_file output_file xmin xmax ymin ymax zmin zmax\n"
         << "Exiting.\n" << endl;
    exit(0);
  }

  try {
    string in_file_name(argv[1]);
    string out_file_name(argv[2]);
    int xmin = atoi(argv[3]);
    int xmax = atoi(argv[4]);
    if (xmax < xmin) throw string("Error xmin > xmax\n");
    int ymin = atoi(argv[5]);
    int ymax = atoi(argv[6]);
    if (ymax < ymin) throw string("Error ymin > ymax\n");
    int zmin = atoi(argv[7]);
    int zmax = atoi(argv[8]);
    if (zmax < zmin) throw string("Error zmin > zmax\n");

    // Read the file
    cerr << "Reading tomogram \""<<in_file_name<<"\"" << endl;
    MrcSimple tomo;
    tomo.Read(in_file_name, false);//You can also use "tomo.Read(cin);" or "cin>>tomo;"
    tomo.PrintStats(cerr);  //Optional (display the tomogram size & format)

    if (xmin < 0) xmin = 0;
    if (ymin < 0) ymin = 0;
    if (zmin < 0) zmin = 0;
    if (xmax > tomo.mrc_header.nvoxels[0]) xmax = tomo.mrc_header.nvoxels[0];
    if (ymax > tomo.mrc_header.nvoxels[1]) ymax = tomo.mrc_header.nvoxels[1];
    if (zmax > tomo.mrc_header.nvoxels[2]) zmax = tomo.mrc_header.nvoxels[2];

    MrcSimple cropped_tomo;
    cropped_tomo.mrc_header = tomo.mrc_header;

    float voxel_size[3];
    voxel_size[0] = cropped_tomo.mrc_header.cellA[0] / cropped_tomo.mrc_header.nvoxels[0];
    voxel_size[1] = cropped_tomo.mrc_header.cellA[1] / cropped_tomo.mrc_header.nvoxels[1];
    voxel_size[2] = cropped_tomo.mrc_header.cellA[2] / cropped_tomo.mrc_header.nvoxels[2];

    // Optional: reduce the cellA[0], cellA[1], cellA[2] size due to cropping

    cropped_tomo.mrc_header.cellA[0] *= (1.0 + xmax-xmin) / cropped_tomo.mrc_header.nvoxels[0];
    cropped_tomo.mrc_header.cellA[1] *= (1.0 + ymax-ymin) / cropped_tomo.mrc_header.nvoxels[1];
    cropped_tomo.mrc_header.cellA[2] *= (1.0 + zmax-zmin) / cropped_tomo.mrc_header.nvoxels[2];
    cropped_tomo.mrc_header.nvoxels[0] = 1 + xmax-xmin;
    cropped_tomo.mrc_header.nvoxels[1] = 1 + ymax-ymin;
    cropped_tomo.mrc_header.nvoxels[2] = 1 + zmax-zmin;
    cropped_tomo.mrc_header.mvoxels[0] = 1 + xmax-xmin;
    cropped_tomo.mrc_header.mvoxels[1] = 1 + ymax-ymin;
    cropped_tomo.mrc_header.mvoxels[2] = 1 + zmax-zmin;
    //           and shift origin[0],origin[1],origin[2] accordingly
    cropped_tomo.mrc_header.origin[0] -= (xmin-1)*voxel_size[0];
    cropped_tomo.mrc_header.origin[1] -= (ymin-1)*voxel_size[1];
    cropped_tomo.mrc_header.origin[2] -= (zmin-1)*voxel_size[2];
    cropped_tomo.Alloc();
    for (int iz = 0; iz <= zmax-zmin; iz++)
      for (int iy = 0; iy <= ymax-ymin; iy++)
	for (int ix = 0; ix <= xmax-xmin; ix++)
	  cropped_tomo.aaafDensity[iz][iy][ix] = 
	    tomo.aaafDensity[iz+zmin][iy+ymin][ix+xmin];

    // Write the file
    cerr << "writing cropped tomogram (in float mode)" << endl;
    cerr << "  new tomogram after cropping:" << endl;
    cropped_tomo.PrintStats(cerr);
    cropped_tomo.Write(out_file_name); //(You can also use cout<<cropped_tomo;)
  }
  catch (string s) {
    cerr << s << endl; // In case of file format error, print message and exit
    exit(1);
  }
}

