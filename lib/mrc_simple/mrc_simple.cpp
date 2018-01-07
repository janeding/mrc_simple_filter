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


#include <cstring>
#include <cassert>
#include <fstream>
#include <iostream>
using namespace std;
#include <alloc3d.h>
#include "mrc_simple.h"


template<class RealNum >
inline RealNum ABS(RealNum x) { return ((x<0.0) ? -x: x); }

void MrcHeader::Read(istream& mrc_file) {
  // Read the first "SIZE_HEADER" bytes into a temporary array "header_data"
  char header_data[MrcHeader::SIZE_HEADER]; //make sure we read the correct number of bytes
  mrc_file.read(header_data, MrcHeader::SIZE_HEADER);

  // copy the bytes we dont use into the remaining_raw_data array
  memcpy(remaining_raw_data,  
         header_data + MrcHeader::SIZE_HEADER_USED, 
         MrcHeader::SIZE_REMAINING_HEADER);
          
  // Convert the binary data into numbers (Ints and floats)

  // 1    NX       number of columns (fastest changing in map)
  // 2    NY       number of rows   
  // 3    NZ       number of sections (slowest changing in map)

  // Commenting out the lines with "mrc_file.read()", although they also work.
  //mrc_file.read((char*)&nvoxels[0], sizeof(Int));
  //mrc_file.read((char*)&nvoxels[1], sizeof(Int));
  //mrc_file.read((char*)&nvoxels[2], sizeof(Int));
  //nvoxels[0] = *((Int*)(pheader_data)+0);
  nvoxels[0] = *(reinterpret_cast<Int*>(header_data)+0);
  nvoxels[1] = *(reinterpret_cast<Int*>(header_data)+1);
  nvoxels[2] = *(reinterpret_cast<Int*>(header_data)+2);


  // 4    MODE     data type :
  //       0        image : signed 8-bit bytes range -128 to 127
  //       1        image : 16-bit halfwords
  //       2        image : 32-bit reals
  //       3        transform : complex 16-bit integers
  //       4        transform : complex 32-bit reals
  //       6        image : unsigned 16-bit range 0 to 65535
  //mrc_file.read((char*)&mode, sizeof(Int));
  mode = *(reinterpret_cast<Int*>(header_data)+3);


  // 5    NXSTART number of first column in map (Default = 0)
  // 6    NYSTART number of first row in map
  // 7    NZSTART number of first section in map
  //mrc_file.read((char*)&nstart[0], sizeof(Int));
  //mrc_file.read((char*)&nstart[1], sizeof(Int));
  //mrc_file.read((char*)&nstart[2], sizeof(Int));
  nstart[0] = *(reinterpret_cast<Int*>(header_data)+4);
  nstart[1] = *(reinterpret_cast<Int*>(header_data)+5);
  nstart[2] = *(reinterpret_cast<Int*>(header_data)+6);


  // 8     MX       number of intervals along X
  // 9     MY       number of intervals along Y
  // 10    MZ       number of intervals along Z
  //mrc_file.read((char*)&mvoxels[0], sizeof(Int));
  //mrc_file.read((char*)&mvoxels[1], sizeof(Int));
  //mrc_file.read((char*)&mvoxels[2], sizeof(Int));
  mvoxels[0] = *(reinterpret_cast<Int*>(header_data)+7);
  mvoxels[1] = *(reinterpret_cast<Int*>(header_data)+8);
  mvoxels[2] = *(reinterpret_cast<Int*>(header_data)+9);


  // 11-13    CELLA    cell dimensions in angstroms
  //mrc_file.read((char*)&(cellA[0]), sizeof(float));
  //mrc_file.read((char*)&(cellA[1]), sizeof(float));
  //mrc_file.read((char*)&(cellA[2]), sizeof(float));
  cellA[0] = *(reinterpret_cast<float*>(header_data)+10);
  cellA[1] = *(reinterpret_cast<float*>(header_data)+11);
  cellA[2] = *(reinterpret_cast<float*>(header_data)+12);


  // 14-16    CELLB    cell angles in degrees
  //mrc_file.read((char*)&(cellB[0]), sizeof(float));
  //mrc_file.read((char*)&(cellB[1]), sizeof(float));
  //mrc_file.read((char*)&(cellB[2]), sizeof(float));
  cellB[0] = *(reinterpret_cast<float*>(header_data)+13);
  cellB[1] = *(reinterpret_cast<float*>(header_data)+14);
  cellB[2] = *(reinterpret_cast<float*>(header_data)+15);


  // 17    MAPC     axis corresp to cols (1,2,3 for X,Y,Z)
  // 18    MAPR     axis corresp to rows (1,2,3 for X,Y,Z)
  // 19    MAPS     axis corresp to sections (1,2,3 for X,Y,Z)
  //mrc_file.read((char*)&(mapCRS[0]), sizeof(Int));
  //mrc_file.read((char*)&(mapCRS[1]), sizeof(Int));
  //mrc_file.read((char*)&(mapCRS[2]), sizeof(Int));
  mapCRS[0] = *(reinterpret_cast<Int*>(header_data)+16);
  mapCRS[1] = *(reinterpret_cast<Int*>(header_data)+17);
  mapCRS[2] = *(reinterpret_cast<Int*>(header_data)+18);


  // 20    DMIN     minimum density value
  // 21    DMAX     maximum density value
  // 22    DMEAN    mean density value
  //mrc_file.read((char*)&(dmin), sizeof(float));
  //mrc_file.read((char*)&(dmax), sizeof(float));
  //mrc_file.read((char*)&(dmean), sizeof(float));
  dmin  = *(reinterpret_cast<float*>(header_data)+19);
  dmax  = *(reinterpret_cast<float*>(header_data)+20);
  dmean = *(reinterpret_cast<float*>(header_data)+21);


  // 23    ISPG     space group number 0 or 1 (default=0)
  //mrc_file.read((char*)&(ispg), sizeof(Int));
  ispg = *(reinterpret_cast<Int*>(header_data)+22);

  // 24    NSYMBT   number of bytes used for symmetry data (0 or 80)
  //mrc_file.read((char*)&(nsymbt), sizeof(Int));
  nsymbt = *(reinterpret_cast<Int*>(header_data)+23);

  // 25-49    EXTRA    extra space used for anything   - 0 by default
  //mrc_file.read(extra_raw_data, 25*MrcHeader::SIZE_PER_FIELD)
  // copy the bytes we dont use into the remaining_raw_data array
  memcpy(extra_raw_data,
         header_data + (24*MrcHeader::SIZE_PER_FIELD),
         25*MrcHeader::SIZE_PER_FIELD);

  // 50-52    ORIGIN   origin in X,Y,Z used for transforms
  //mrc_file.read((char*)&(origin[0]), sizeof(float));
  //mrc_file.read((char*)&(origin[1]), sizeof(float));
  //mrc_file.read((char*)&(origin[2]), sizeof(float));
  origin[0] = *(reinterpret_cast<float*>(header_data)+49);
  origin[1] = *(reinterpret_cast<float*>(header_data)+50);
  origin[2] = *(reinterpret_cast<float*>(header_data)+51);

  if (mode == MrcHeader::MRC_MODE_BYTE) {
    // Try to determine if 8-bit integers are signed or unsigned.
    // Apparently, there are no general guidelines indicating whether 8-bit 
    // integers are signed or unsigned in MRC files.  
    // If either of these statements are true, bytes are signed.
    //   1) check the filename to see if it ends in ".rec" instead of ".mrc"
    //     http://www.cgl.ucsf.edu/pipermail/chimera-users/2010-June/005245.html
    //   2) check for a bit in the "imodFlags" section
    //      of the header.  (That's what we do here.)
    //     http://bio3d.colorado.edu/imod/doc/mrc_format.txt
    // Otherwise, the are (by default) assumbe to be unsigned.
    //
    //152 230  4   int     imodStamp   1146047817 indicates that file was created by IMOD or 
    //                                 other software that uses bit flags in the following field
    //156 234  4   int     imodFlags   Bit flags:
    //                                 1 = bytes are stored as signed
    //                                 2 = pixel spacing was set from size in extended header 
    //                                 4 = origin is stored with sign inverted from definition
    //                                     below
    Int imodStamp = *(reinterpret_cast<Int*>(header_data)+38); //(byte 152-155)
    if (imodStamp == 1146047817) { 
      // Then the file was created by IMOD or other software using "bit flags"
      // See: http://bio3d.colorado.edu/imod/doc/mrc_format.txt
      Int imodFlags = *(reinterpret_cast<Int*>(header_data)+39);//(byte 156-159)
      use_signed_bytes = (imodFlags & 1);
    }
  }

} //MrcHeader::Read(istream& mrc_file)




void MrcHeader::Write(ostream& mrc_file) const {

  // 1    NX       number of columns (fastest changing in map)
  // 2    NY       number of rows   
  // 3    NZ       number of sections (slowest changing in map)

  // Commenting out the lines with "mrc_file.read()", although they also work.
  mrc_file.write((char*)&nvoxels[0], sizeof(Int));
  mrc_file.write((char*)&nvoxels[1], sizeof(Int));
  mrc_file.write((char*)&nvoxels[2], sizeof(Int));

  // 4    MODE     data type :
  //       0        image : signed 8-bit bytes range -128 to 127
  //       1        image : 16-bit halfwords
  //       2        image : 32-bit reals
  //       3        transform : complex 16-bit integers
  //       4        transform : complex 32-bit reals
  //       6        image : unsigned 16-bit range 0 to 65535
  mrc_file.write((char*)&mode, sizeof(Int));


  // 5    NXSTART number of first column in map (Default = 0)
  // 6    NYSTART number of first row in map
  // 7    NZSTART number of first section in map
  mrc_file.write((char*)&nstart[0], sizeof(Int));
  mrc_file.write((char*)&nstart[1], sizeof(Int));
  mrc_file.write((char*)&nstart[2], sizeof(Int));


  // 8    MX       number of intervals along X
  // 9    MY       number of intervals along Y
  // 10    MZ       number of intervals along Z
  mrc_file.write((char*)&mvoxels[0], sizeof(Int));
  mrc_file.write((char*)&mvoxels[1], sizeof(Int));
  mrc_file.write((char*)&mvoxels[2], sizeof(Int));


  // 11-13    CELLA    cell dimensions in angstroms
  mrc_file.write((char*)&(cellA[0]), sizeof(float));
  mrc_file.write((char*)&(cellA[1]), sizeof(float));
  mrc_file.write((char*)&(cellA[2]), sizeof(float));


  // 14-16    CELLB    cell angles in degrees
  mrc_file.write((char*)&(cellB[0]), sizeof(float));
  mrc_file.write((char*)&(cellB[1]), sizeof(float));
  mrc_file.write((char*)&(cellB[2]), sizeof(float));


  // 17    MAPC     axis corresp to cols (1,2,3 for X,Y,Z)
  // 18    MAPR     axis corresp to rows (1,2,3 for X,Y,Z)
  // 19    MAPS     axis corresp to sections (1,2,3 for X,Y,Z)
  mrc_file.write((char*)&(mapCRS[0]), sizeof(Int));
  mrc_file.write((char*)&(mapCRS[1]), sizeof(Int));
  mrc_file.write((char*)&(mapCRS[2]), sizeof(Int));

  // 20    DMIN     minimum density value
  // 21    DMAX     maximum density value
  // 22    DMEAN    mean density value
  mrc_file.write((char*)&(dmin), sizeof(float));
  mrc_file.write((char*)&(dmax), sizeof(float));
  mrc_file.write((char*)&(dmean), sizeof(float));


  // 23    ISPG     space group number 0 or 1 (default=0)
  mrc_file.write((char*)&(ispg), sizeof(Int));

  // 24    NSYMBT   number of bytes used for symmetry data (0 or 80)
  mrc_file.write((char*)&(nsymbt), sizeof(Int));

  // 25-49    EXTRA    extra space used for anything   - 0 by default
  mrc_file.write(extra_raw_data, (1+(49-25))*sizeof(Int));

  // 50-52    ORIGIN   origin in X,Y,Z used for transforms
  mrc_file.write((char*)&(origin[0]), sizeof(float));
  mrc_file.write((char*)&(origin[1]), sizeof(float));
  mrc_file.write((char*)&(origin[2]), sizeof(float));

  // I don't care about the remaining junk in the header
  mrc_file.write(remaining_raw_data,
                 MrcHeader::SIZE_REMAINING_HEADER);

} //MrcHeader::Write(ostream& mrc_file) const





void MrcHeader::PrintStats(ostream& out) {
  out << "tomogram number of voxels ("<< nvoxels[0] << ", " 
                                      << nvoxels[1] << ", "
                                      << nvoxels[2] << ")" 
                                      << endl;
  //out << "tomogram physical dimensions ("<< cellA[0] << "," << cellA[1] << "," << cellA[2] << ")" << endl;
  out << "tomogram voxel size in nm ("
      << 0.1*cellA[0]/nvoxels[0] << ", " 
      << 0.1*cellA[1]/nvoxels[1] << ", " 
      << 0.1*cellA[2]/nvoxels[2] << ")" << endl;
  out << "tomogram table axis order (" << mapCRS[0] << ", " 
                                       << mapCRS[1] << ", "
                                       << mapCRS[2] << ")"
                                       << endl;
  out << "tomogram mode " << mode << endl;
  out << "tomogram minimum, maximum density: " << dmin << " "
                                               << dmax<<endl;
  out << "tomogram origin ("<< origin[0] << ", " 
                            << origin[1] << ", " 
                            << origin[2] << ")" << endl;
}





void MrcSimple::Alloc() {
  Alloc3D(mrc_header.nvoxels,
          &afDensity,
          &aaafDensity);
}


void MrcSimple::Dealloc() {
  Dealloc3D(mrc_header.nvoxels,
            &afDensity,
            &aaafDensity);
}





void MrcSimple::Read(istream& mrc_file,
                     bool rescale,
		     float ***aaafMask) {
  mrc_header.Read(mrc_file);
  if ((mrc_header.mapCRS[0] != 1) ||
      (mrc_header.mapCRS[1] != 2) ||
      (mrc_header.mapCRS[2] != 3))
    throw "Error: This program currently only supports .MRC/.REC files in row-major\n"
      "       format, and whose mapC, mapR, mapS numbers (in the file header)\n"
      "       are 1,2,3, respectively.  Please use an alternate program to convert\n"
      "       to this format.  For more information on MRC file format, see:\n"
      "       http://www2.mrc-lmb.cam.ac.uk/research/locally-developed-software/image-processing-software/#image"
      "       http://bio3d.colorado.edu/imod/doc/mrc_format.txt\n"
      "       http://ami.scripps.edu/software/mrctools/mrc_specification.php\n";
  ReadArray(mrc_file);
  if (rescale)
    Rescale01(aaafMask);
}



void MrcSimple::Read(string in_file_name,
                     bool rescale,
		     float ***aaafMask) {
  Int len_in_file_name = in_file_name.size();
  fstream mrc_file;
  mrc_file.open(in_file_name.c_str(), ios::binary | ios::in);
  if (! mrc_file) 
    throw "Error: unable to open \"" + in_file_name + "\" for reading.\n";
  // Try to infer signed-vs-unsigned integers from the file name:
  //http://www.cgl.ucsf.edu/pipermail/chimera-users/2010-June/005245.html
  if ((len_in_file_name > 4)
      && 
      (in_file_name.substr(len_in_file_name-4, len_in_file_name) == ".rec")) {
    mrc_header.use_signed_bytes = false;
  }
  Read(mrc_file, rescale, aaafMask);
  mrc_file.close();
}


void MrcSimple::ReadArray(istream& mrc_file) {
  Dealloc(); //free up any space you may have allocated earlier
  Alloc();   //allocate space for the array
  //mrc_file.read((char*)afDensity, sizeof(float)*num_voxels);
  for(Int iz=0; iz<mrc_header.nvoxels[2]; iz++) {
    //aaafDensity[iz] = new float* [mrc_header.nvoxels[1]];
    for(Int iy=0; iy<mrc_header.nvoxels[1]; iy++) {
      //aaafDensity[iz][iy] = &(afDensity[iz*mrc_header.nvoxels[0]*mrc_header.nvoxels[1] + 
      //                                    iy*mrc_header.nvoxels[0]]);
      for(Int ix=0; ix<mrc_header.nvoxels[0]; ix++) {

        switch (mrc_header.mode) {

        case MrcHeader::MRC_MODE_BYTE:
          if (mrc_header.use_signed_bytes)
          {
            int8_t entry;
            mrc_file.read((char*)&entry, sizeof(char));
            aaafDensity[iz][iy][ix] = static_cast<float>(entry);
          }
          else
          {
            uint8_t entry;
            mrc_file.read((char*)&entry, sizeof(char));
            aaafDensity[iz][iy][ix] = static_cast<float>(entry);
          }
          break;

        case MrcHeader::MRC_MODE_SHORT:
          {
            int16_t entry;
            mrc_file.read((char*)&entry, sizeof(int16_t));
            aaafDensity[iz][iy][ix] = static_cast<float>(entry);
          }
          break;

        case MrcHeader::MRC_MODE_USHORT:
          {
            uint16_t entry;
            mrc_file.read((char*)&entry, sizeof(uint16_t));
            aaafDensity[iz][iy][ix] = static_cast<float>(entry);
          }
          break;

        case MrcHeader::MRC_MODE_FLOAT:
          {
            float entry;
            mrc_file.read((char*)&entry, sizeof(float));
            //aaafDensity[iz][iy][ix] = static_cast<float>(entry);
            aaafDensity[iz][iy][ix] = entry;  //(already type float)
          }
          break;
        default:
          throw string("UNSUPPORTED MODE in MRC file (unsupported MRC format)");
          exit(-1);
          break;
        } // switch (mrc_header.mode)

        // Debugging check:
        //if ((abs(ix-100) <= 2) && (iy==100) && (iz==100))
        //  cerr << "aaafDensity["<<iz<<"]["<<iy<<"]["<<ix<<"] = "
        //       << aaafDensity[iz][iy][ix] << endl;

      } //for(Int ix=0; ix<mrc_header.nvoxels[0]; ix++) {
    } //for(Int iy=0; iy<mrc_header.nvoxels[1]; iy++) {
  } //for(Int iz=0; iz<mrc_header.nvoxels[2]; iz++) {
} //MrcSimple::ReadArray(istream& mrc_file)



void MrcSimple::Write(string out_file_name) {
  fstream mrc_file;
  mrc_file.open(out_file_name.c_str(), ios::binary | ios::out);
  if (! mrc_file) 
    throw "Error: unable to open \"" + out_file_name + "\" for writing.\n";
  Write(mrc_file);  // You can also use "mrc_file << tomo;"
  mrc_file.close();
}


void MrcSimple::Write(ostream& mrc_file) {
  // First, write the MRC file header:
  //mrc_header.Write(mrc_file); <-- OOPS, that does not work.  Commenting out.
  // Note that any changes made to the tomogram density data will effect the
  // header. So we make a copy of the original header and modify it accordingly:
  // As of 2015-4-16, regardless of the original header, I save the list
  // of numbers in "float" format.  This means I must change mrc.mode.
  MrcHeader new_mrc_header = mrc_header;
  new_mrc_header.mode = MrcHeader::MRC_MODE_FLOAT;

  FindMinMaxMean(); // calculates the "dmin", "dmax", and "dmean" member values
  new_mrc_header.Write(mrc_file);

  // Finally, write the array of densities:
  WriteArray(mrc_file);

} //MrcSimple::Write(ostream& mrc_file) const






void MrcSimple::WriteArray(ostream& mrc_file) const {
  // Write the MRC contents:
  for(Int iz=0; iz<mrc_header.nvoxels[2]; iz++)
    for(Int iy=0; iy<mrc_header.nvoxels[1]; iy++)
      for(Int ix=0; ix<mrc_header.nvoxels[0]; ix++)
        mrc_file.write((char*)&(aaafDensity[iz][iy][ix]), sizeof(float));
}



void MrcSimple::FindMinMaxMean(float ***aaafMask) {
  double density_total = 0.0;
  double dmin = 0.0;  //impossible initial value
  double dmax = -1.0; //impossible initial value
  long long nvoxels_kept = 0;
  for(Int iz=0; iz<mrc_header.nvoxels[2]; iz++) {
    for(Int iy=0; iy<mrc_header.nvoxels[1]; iy++) {
      for(Int ix=0; ix<mrc_header.nvoxels[0]; ix++) {
	if (aaafMask && (aaafMask[iz][iy][ix] == 0))
	  continue; //if a mask is supplied, ignore voxels when mask=0
        density_total += aaafDensity[iz][iy][ix];
        if (dmin > dmax) {
          dmin = aaafDensity[iz][iy][ix];
          dmax = aaafDensity[iz][iy][ix];
        }
        else {
          if (aaafDensity[iz][iy][ix] > dmax)
            dmax = aaafDensity[iz][iy][ix];
          if (aaafDensity[iz][iy][ix] < dmin)
            dmin = aaafDensity[iz][iy][ix];
        }
	nvoxels_kept++;
      }
    }
  }
  mrc_header.dmin = dmin;
  mrc_header.dmax = dmax;
  mrc_header.dmean = density_total  /  nvoxels_kept;
}


void MrcSimple::Rescale01(float ***aaafMask)
{
  FindMinMaxMean(aaafMask); // set mrc_header.dmin,dmax,dmean
  float dmin = mrc_header.dmin;
  float dmax = mrc_header.dmax;
  for(Int iz=0; iz<mrc_header.nvoxels[2]; iz++) {
    for(Int iy=0; iy<mrc_header.nvoxels[1]; iy++) {
      for(Int ix=0; ix<mrc_header.nvoxels[0]; ix++) {
	if (aaafMask && (aaafMask[iz][iy][ix] == 0))
	  continue; //if a mask is supplied, ignore voxels when mask=0
        aaafDensity[iz][iy][ix] = 
          (aaafDensity[iz][iy][ix] - dmin) / (dmax - dmin);
      }
    }
  }
  mrc_header.dmean = (mrc_header.dmean - dmin) / (dmax - dmin);
  mrc_header.dmin = 0.0;
  mrc_header.dmax = 1.0;
}



void MrcSimple::Invert(float ***aaafMask)
{
  double sum = 0.0;
  long n = 0;
  for(Int iz=0; iz<mrc_header.nvoxels[2]; iz++) {
    for(Int iy=0; iy<mrc_header.nvoxels[1]; iy++) {
      for(Int ix=0; ix<mrc_header.nvoxels[0]; ix++) {
	if (aaafMask && (aaafMask[iz][iy][ix] == 0))
	  continue; //if a mask is supplied, ignore voxels when mask=0
        sum += aaafDensity[iz][iy][ix];
        n += 1;
      }
    }
  }
  double ave = sum / n;
  float dmin = ave;
  float dmax = ave;
  for(Int iz=0; iz<mrc_header.nvoxels[2]; iz++) {
    for(Int iy=0; iy<mrc_header.nvoxels[1]; iy++) {
      for(Int ix=0; ix<mrc_header.nvoxels[0]; ix++) {
	if (aaafMask && (aaafMask[iz][iy][ix] == 0))
	  continue; //if a mask is supplied, ignore voxels when mask=0
        aaafDensity[iz][iy][ix] = 2.0*ave - aaafDensity[iz][iy][ix];
        if (aaafDensity[iz][iy][ix] < dmin)
          dmin = aaafDensity[iz][iy][ix];
        if (dmax < aaafDensity[iz][iy][ix])
          dmax = aaafDensity[iz][iy][ix];
      }
    }
  }
  //mrc_header.dmean remains the same, and should equal "ave"
  //(Note: I could also use FindMinMaxMean(aaafMask))
  assert(ABS(mrc_header.dmean - ave) < 1.0e-3*ABS(ave));
}

