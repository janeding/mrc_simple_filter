convolve_mrc
===========

**convolve_mrc** applies a filter to a tomogram in the X,Y,Z directions
and saves the result as a new .mrc/.rec file.

This program supports masks ans thresholding using
the "*-mask*" and "*-thresh*" (and "*-thresh2*" and "*-thresh4*") arguments.

### Distance measurements:

In this program all of the distances provided by the user are assumed to be in *physical units*.  (Not *voxels*.  This effects all of the gaussian "sigma" parameters, specifically all of the "-width-xy" and "-width-z" and "-window" arguments.  To express these distances in *voxels* instead of *nm*, set the voxel_width to 1.0 using "-w 1.0".)


## Usage Example:

```
convolve_mrc \
   -in Bdellovibrio_1K.rec \
   -out Bdellovibrio_1K_ribosome_peaks.mrc \
   -w 1.92 \
   -dog 12.0 12.0 24.0 15.0 15.0 30.0 \
   -exponents 2 2 \
   -cutoff 2.0 \
   -mask Bdellovibrio_1K_mask_water=0_periplasm=1_cytoplasm=2.mrc \
   -mask-select 0 -mask-out 1.0
```

## Arguments:

## Input and Output files

The user must specify the name of the tomogram they wish to process using the
"-in" argument:
```
   -in SOURCE_FILE.mrc
```
(Note: files may also end in ".rec")

Users must specify the name of the new tomogram
(created by applying the filter to the original tomogram)
using the "-out" argument:
```
   -out DESTINATION_FILE.mrc
```

## Voxel Width

```
   -w voxelwidth
```
This optional argument allows you to manually specify the width of each voxel
in *physical units* (as opposed to, *number-of-voxels*).
This is necessary because all distance parameters entered by the user must be expressed in physical units, not voxels.
If not specified,
*by default*
the physical width of each voxel will be read from the MRC file.

*(Using the "-w" argument will override the voxel widths specified in the header of the MRC file.
Note: In most of the Jensen lab tomograms before 2018,
 voxel widths in the tomogram files are expressed units of Angstroms
 not nm.  So if you do not use the -w argument, you should keep this in
 mind when you specify your filter width parameters.
 Also note that there is typically a ~10% difference between the voxel
 width stored in the header of an MRC file, and the voxel width stored at the
 corresponding entry in the tomography database.  The later is more accurate.)*

## Selecting a Filter:

The user can select the type of filter they want to use
using command line arguments.
The "-gauss" filter uses a (low-pass) Gaussian filter. (default)
The "-dog" filter uses a (band-pass) Difference-Of-Gaussians filter
which can be used for blob detection as well as low-frequency removal.
The "-dogxy" filter uses a Difference-Of-Gaussians filter in the XY direction
coupled with an ordinary Gaussian filter in the Z direction.
In these examples, the Gaussians can be customized,
although it may slow down the filtering process significantly.

## -gauss
The -gauss argument must be followed by 3 numbers:
```
   -gauss  s_x  s_y  s_z
```
If the "-gauss" filter is selected, the
original image is convolved with the following function:
```
   h(x,y,z) = A*exp(-0.5 * r^2)
    where r = sqrt((x/s_x)^2 + (y/s_y)^2 + (z/s_z)^2))
```
The width of the Gaussian (ie, the s_x,s_y,s_z arguments) should be specified in units of physical distance, not in voxels.
By default the domain of the filter is extended in each direction to a distance
of 3.0 * MAX(a,b).  (IE, thrice the width of the filter in that direction. This can be overridden using the "-cutoff" argument.  See below.)
(The A coefficient will determined automatically by normalization, ie., so that the discrete sum of h(x,y,z) over x,y,z is 1.)

If the optional "-exponent n" argument is supplied, then
the original image is instead convolved with the following function:
```
h(x,y,z) = A*exp(-r^n)
 where r = sqrt((x/s_x)^2 + (y/s_y)^2 + (z/s_z)^2))
```
This will slow down the filter considerably.
(Because, in this case the filter is no longer a separable function of x,y,z
 and the full 3-D convolution must be performed.)


## -dog
The -dog argument must be followed by 6 numbers:

```
  -dog  a_x  a_y  a_z  b_x  b_y  b_z
```
If the "-dog" filter is selected, the
original image is convolved with the following function:
```
   h(x,y,z) = A*exp(-0.5 * r_a^2) - B*exp(-0.5 * r_b^2)
  where r_a = sqrt((x/a_x)^2 + (y/a_y)^2 + (z/a_z)^2))
    and r_b = sqrt((x/b_x)^2 + (y/b_y)^2 + (z/b_z)^2))
```
The width of the Gaussian (the a_x,a_y,a_z,b_x,b_y,b_z arguments) should be specified in units of physical distance, not in voxels.
The A and B coefficients will be automatically chosen so that the discrete sum of h(x,y,z) over x,y,z is 0, and the peak height is 1 (A-B=1).
By default the domain of the filter is extended in each direction to a distance
of 3.0 * MAX(a,b).  (IE, thrice the width of the filter in that direction. This can be overridden using the "-cutoff" argument.  See below.)

If the optional "-exponents m n" argument is supplied, then
the original image is instead convolved with the following function:
```
h(x,y,z) = A*exp(-r_a^m) - B*exp(-r_b^n)
 where r_a = sqrt((x/a_x)^2 + (y/a_y)^2 + (z/a_z)^2))
   and r_b = sqrt((x/b_x)^2 + (y/b_y)^2 + (z/b_z)^2))
```
This will slow down the filter considerably.
This crude generalization of the "gaussian" function gives us an ad-hoc way
to alter shape of the filter.


## -dogxy
The -dogxy argument must be followed by 3 numbers:
```
  -dog  a  b  c
```
If the "-dogxy" filter is selected,
the original image is convolved with the following function:
```
   h(x,y,z) = h_xy(x,y) * h_z(z)
```
 In the XY plane, the filter used is:
```
   h_xy(x,y) = A*exp(-(|r|/a)^m) - B*exp(-(|r|/b)^n)
           r = sqrt(x^2 + y^2)
 and A,B are chosen so that discrete sum over x,y of h_xy(x,y) is 0, and A-B=1

```
 The "m" and "n" parameters are exponents.
 They are both set to 2 by default, however you can override them
 using the "-exponents m n" argument.

 Along the Z direction, the filter used is a simple Gaussian:
```
   h_z(r) = C*exp(-0.5*(z/c)^2)
```

## Additional Arguments:

## Filter Size
```
   -cutoff ratio
```
This specifies the number of voxels in the filter which will be convolved
with the image.  
It is expressed in units of the "width" parameters
for the filter you have selected.
For example, if you use "-cutoff 2.0" with the "-gauss" (Gaussian) filter,
then the filter window will extend outward away from the filter center
by a distance of 2 sigma in the respective x,y,z directions.
(ie., 2.0*a_x, 2.0*a_y, 2.0*a_z, where "a_x", "a_y", "a_z"
 are the Gaussian widths in the x,y,z directions, respectively).
If unspecified, this parameter is set to 2.0.

```
   -window Wx Wy Wz
```
If you prefer to specify the size of the filter window manually, you can
use the "-window" argument.
This argument specifies the size of the 3-D filter box in x,y, directions
(in physical units, not voxels).
This overrides the "-cutoff ratio" argument.

*(Note: Either way, in the worst case, the product of the 3 numbers, Wx*Wy*Wz, is proportional to the running time of the filter.
However when using the "-gauss" or "-dog" filters with default settings,
the running time is proportional to the sum of these numbers, Wx+Wy+Wz.)*


```
   -ang-to-nm
```
When reading MRC files, this argument will
multiply the voxel width
specified in the MRC file by 0.1.
This allows you to specify your distance parameters in units of nm
instead of Angstroms.

This argument has no effect if the "-w" argument is used.



## Masking

The optional "-mask", "-mask-select", and "-mask-out" arguments allow you to
ignore certain voxels from the source image (tomogram).

Using "masks" allows the user to perform the filtering considering only a
subset of voxels in the tomogram (ie, the "mask").

```
   -mask  file.mrc
```
The argument following the
"-mask" argument should be the name of a tomogram file (MRC format) of the
same size as the input tomogram.  By default only voxels belonging to the
tomogram with non-zero voxel values in the "mask" tomogram will be considered.
(This means if, during the filtering process, part of the filter
falls outside the mask, then those voxels are excluded from the sum, and the
weights from the remaining voxels will be increased (normalized) accordingly.)
```
   -mask-out  brightness_value
```
If the "-mask-out" argument is specified, then voxels outside the mask will
be assigned to the number following this argument.  (Otherwise they are
assigned to 0 by default.)
```
   -mask-select  brightness_value
```
If the "-mask-select" argument is specified, then instead of considering all
voxels with non-zero values from the "mask" image/tomogram, only voxels whose mask
value equal the number following this argument will belong to the mask.  (WARNING:  Do not supply a number outside the range from 0 to 1.  Voxels from 8-bit or 16-bit integer images/tomograms are rescaled from 0 to 1 before performing this operation.  So, in practice, this means that only values of "0" or "1" can be used with abit or 16bit images/tomograms.)

## Rescaling and Thresholding

The optional rescaling and thresholding arguments provide a way to make sure that the brightnesses in the resulting image lie between 0 and 1.
*(Again, Voxels from 8-bit or 16-bit integer images are interpreted as floating point numbers between 0 and 1 by dividing them by 255 and 65535, respectively.)*
The result of the filtering operations above may produce images/tomograms with
brightness values which may like outside the range of values accepted
by other software.
The following arguments provide different ways to rescale the brightness
of each voxel between 0 and 1 according to its brightness value.

```
   -thresh  thresh01
```

If the "-thresh" argument is passed, it must be followed by a number ("thresh01").  *After* the filter is applied, voxels with brightnesses below this number will be replaced with a voxel of brightness 0, and voxels above this number will be replaced with 1.

```
 output
 brightness
 (a.k.a. "density")
  /|\                      _____________________________\
 1 |                      |                             /
   |                      |                 
   |                      |            
 0 |______________________|  ___________________________\ input
                     threshold                          / brightness
                                                          (a.k.a. "density")
```

If the "-thresh2" argument is passed, then it must be followed by 2 numbers:

```
   -thresh2  thresh_01_a  thresh_01_b
```

In this case, *after* the filter is applied,
the resulting voxel brightnesses will be scaled
between 0 and 1 according to the following function:

### if thresh_01_a < thresh_01_b
```
 output
 brightness
 (a.k.a. "density")
  /|\                              _________________\
 1 |                           _.-'                 /
   |                       _,-'                 
   |                   _,-'            
 0 |________________,-'                     ________\ input
                 thresh         thresh              / brightness
                  01_a           01_b                 (a.k.a. "density")
```

***Or***, if ***thresh_01_b < thresh_01_a***, then the output is inverted:

```
 output
 brightness
 (a.k.a. "density")
  /|\
   |________________
 1 |                `-._
   |                    `-._
   |                        '-._            
 0 |                            `-.___________________\ input
                thresh         thresh                 / brightness
                 01_b           01_a                    (a.k.a. "density")
```


If the "-thresh4" argument is passed, then it must be followed by 4 numbers:

```
   -thresh4  thresh01_a  thresh01_b  thresh10_a  thresh10_b
```

In this case, the resulting image voxels will be scaled
between 0 and 1 according to the following function:

```
 output
 brightness
 (a.k.a. "density")
  /|\
 1 |                 ________________                
   |             _,-'                `-._
   |         _,-'                        `-._
 0 |______,-'                                `-._______\ input
        thresh    thresh          thresh    thresh     / brightness
         01_a      01_b            10_a      10_b        ("density")
```
***Or***, if the user selects ***thresh_10_b < thresh_01_a***, then the output is inverted:
```
  /|\                                                   
 1 |_____                                       _______\ input
   |     `-._                               _.-'       / brightness
   |         `-._                       _,-'             ("density")
 0 |             `-._________________,-'
        thresh    thresh          thresh    thresh
         10_a      10_b            01_a      01_b
```
