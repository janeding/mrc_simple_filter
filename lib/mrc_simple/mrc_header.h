#ifndef MRC_HEADER_H
#define MRC_HEADER_H

#include <cstdint>
using namespace std;


// "MrcSimple": a class to read and write tomographic data.
//              (stored in .MRC and .REC files)
//
// Tomographic data is internally stored as arrays of floats, however
// this program can read MRC/REC files using several other numeric formats 
// as well (including signed and unsigned 8-bit and 16-bit integers).
// Note: As of 2015-4-20, file data is assumed to be in row-major format.
//       http://en.wikipedia.org/wiki/Row-major_order
//       (IE The "mapC", "mapR", and "mapS" (mapCRS[]) header data is ignored.)
// Note: As of 2015-4-20, this software does not attempt to detect 
//       or convert big-endian or little-endian numbers.
//       To be on the safe side, compile this software on hardware which is
//       compatible with the hardware which was used to write the MRC file.
//       (IE intel CPUs, or ARM CPUs using the same endian settings.)
//
// Documentation on the MRC file format is available here:
//http://bio3d.colorado.edu/imod/doc/mrc_format.txt
//http://www2.mrc-lmb.cam.ac.uk/research/locally-developed-software/image-processing-software/#image
//http://ami.scripps.edu/software/mrctools/mrc_specification.php
//http://bio3d.colorado.edu/imod/betaDoc/libhelp/mrcfiles.html#TOP




//typedef int Int;
typedef int32_t Int;  //An integer type gauranteed to be exactly 32 bits wide.
                      //(All integers in the header are 32 bits wide)
                      //(As of 2015-3, if you use gcc/g++, you must compile 
                      // using the -std=c++11 flag to support "int32_t")






class MrcHeader {

public:

  Int nvoxels[3];  // How many voxels are in the tomogram? (in x,y,z direction)
                   // (total number = nvoxels[0]*nvoxels[1]*nvoxels[2])

  Int mode; // Format of the numbers representing the density at each voxel 
            // location in the array (as stored in an MRC/REC file):
  const static Int  MRC_MODE_BYTE          = 0;
  const static Int  MRC_MODE_SHORT         = 1;
  const static Int  MRC_MODE_FLOAT         = 2;
  const static Int  MRC_MODE_COMPLEX_SHORT = 3;  // (IGNORED as of 2015-4-16)
  const static Int  MRC_MODE_COMPLEX_FLOAT = 4;  // (IGNORED as of 2015-4-16)
  const static Int  MRC_MODE_USHORT        = 6;
  const static Int  MRC_MODE_RGB           = 16; // (IGNORED as of 2015-4-16)
            // Note: This is only useful when reading or writing MRC files.
            // (Internally, densities are represented using arrays of floats.)


  // ----- (The "nstart" and "mvoxels" entries are not useful for us, but -----
  // -----  they appear early in the header file so I include them here.) -----
  Int nstart[3];   // starting-point for sub-images (IGNORED)             -----
  Int mvoxels[3];  // not sure what this is (IGNORED)                     -----
  // --------------------------------------------------------------------------


  float cellA[3]; // physical size (in Angstroms) of the rectangular box 
                  // containing the tomogram.  Useful to calculate voxel size:
                  // (voxel-width in x direction = cellA[0]/nvoxels[0]
                  //  voxel-width in y direction = cellA[1]/nvoxels[1]
                  //  voxel-width in z direction = cellA[2]/nvoxels[2])

  void Read(istream& mrc_file); //Read an .MRC/.REC file header
  void Write(ostream& mrc_file) const; //Write an .MRC/.REC file header

  // constructor sets the tomogram size to 0,0,0
  MrcHeader() {
    nvoxels[0] = 0; 
    nvoxels[1] = 0; nvoxels[2] = 0; mode = -1; 
    use_signed_bytes = true;
  }






  // ----------------------------------------------------------
  // ------------------- Ignore the rest ----------------------
  // ----------------------------------------------------------





  // (The remaining entries in a header file are not currently useful to me.)

  float cellB[3];  // unit cell angles (IGNORED, assume 90 degrees)

  // As of 2015-4-16, mapCRS numbers currently IGNORED (assumed to be 1,2,3)
  Int mapCRS[3];
  // mapCRS[0] = which axis (x,y,z) corresponds to the "columns" in the array
  // mapCRS[1] = which axis (x,y,z) distinguishes different "rows" in the array
  // mapCRS[2] = which axis (x,y,z) distinguishes the sections (planes) in the array

  float dmin;           // Minimum-possible numeric value of entry in the array
  float dmax;           // Maximum-possible numeric value of entry in the array
  float dmean;          // Average value of numbers in the array (IGNORED)
  Int ispg;             // no idea what this is (IGNORED)
  Int nsymbt;           // no idea what this is (IGNORED)

  static const Int SIZE_HEADER = 1024;
  static const Int SIZE_PER_FIELD = 4;  //Each Int or float requires 4 bytes
  static const Int NUM_USED_FIELDS = 52; //Number of entries read from header
  static const Int SIZE_HEADER_USED = NUM_USED_FIELDS*SIZE_PER_FIELD;//in bytes
  static const Int SIZE_REMAINING_HEADER = SIZE_HEADER - SIZE_HEADER_USED;

  char extra_raw_data[25*SIZE_PER_FIELD]; // (IGNORED)

  float origin[3];      // I'm guessing this is only relevant to reconstructions
                        // I'm guessing this might be the location of a point
                        // within the pivot axis (axes?) around which the sample
                        // was rotated when generating the original tilt series.

  char remaining_raw_data[SIZE_REMAINING_HEADER];

  bool use_signed_bytes;// signed bytes? (only relevant when mode=MRC_MODE_BYTE)
  void PrintStats(ostream &out);//Print information about tomogram size & format
};




// Although not necessarily recommended, you can use << or >> to read/write
inline istream& operator>>(istream& mrc_file, 
                           MrcHeader& mrc_header)
{
  mrc_header.Read(mrc_file);
  return mrc_file;
}


inline ostream& operator<<(ostream& mrc_file, 
                           const MrcHeader& mrc_header)
{
  mrc_header.Write(mrc_file);
  return mrc_file;
}




#endif //#ifndef MRC_HEADER_H


