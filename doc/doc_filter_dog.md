filter_dog
===========

**filter_dog** applies a "difference-of-gaussians" filter to the tomogram
in the XY directions, and an ordinary gaussian filter in the Z direction.
It saves the result as a new .mrc/.rec file.
This program was intended to be used to detect ribosomes and other blobs in electron tomography images.

As with "*filter_gauss*", this program supports masks ans thresholding using
the "*-mask*" and "*-thresh*" (and "*-thresh2*" and "*-thresh4*") arguments.
Detailed documentation for this program is located in the "*doc/*" subdirectory.

### Distance measurements:

In this program all of the distances provided by the user are assumed to be in *nm*.  (Not *voxels*.  This effects all of the gaussian "sigma" parameters, specifically all of the "-width-xy" and "-width-z" and "-window" arguments.  To express these distances in *voxels* instead of *nm*, set the voxel_width to 1.0 using "-w 1.0".)


## Usage Example:

```
./filter_dog \
   -in Bdellovibrio_1K.rec \
   -width-xy 12.0 15.0 -width-z 3.0 \
   -cutoff-xy 2.0 -cutoff-z 1.5 \
   -exponents 4 4 \
   -w 1.92 \
   -mask Bdellovibrio_1K_mask_water=0_periplasm=1_cytoplasm=2.mrc \
   -mask-select 0 -mask-out 1.0 \
   -out Bdellovibrio_1K_ribosome_peaks.mrc
```

## Filter Used:
The original image brightness is convolved with the following function:
```
   h(x,y,z) = h_xy(r) * h_z(z)
```
 In the XY plane, the filter used is:
```
   h_xy(r) = A*exp(-(|r|/s)^m) - B*exp(-(|r|/t)^n)
```
 Here, "r" refers to disance from the Z axis
 "s" and "t" are widths in nm (specified by the arguments following "-width-xy")
 "m" and "n" are exponents (specified by the arguments following "-exponents")
 and "A" and "B" are chosen automatically so that h_xy(r=0) = A+B = 1, and
 the sum of h_xy(r) for all pixels in the filter range is zero.
 This crude generalization of the "gaussian" function gives us an ad-hoc way to alter the shape of the filter.

 Along the Z direction, the filter used is a simple Gaussian:
```
   h_z(r) = C*exp(-0.5*(z/u)^2)
```
 Where "u" is specified by the argument following "-width-z"
 and "C" is chosen so that the (discrete) sum of h_z(z) values along z is 1


## Arguments:

## Voxel Width
```
   -w voxelwidth
```
This argument allows you to specify the width of each voxel (in *nm*).
It is assumed that the width of each voxel in the x,y,z directions is the same.

## Filter Size
```
   -window Wx Wy Wz
```
This argument specifies the size of the 3-D filter box (in nm, not voxels).
The product of the 3 numbers (Wx*Wy*Wz) is proportional to the running
time of the filter.
```
   -cutoff-xy  distanceXY
```
This specifies the width of the filter used in the XY direction
(in units of the "t" parameter from the formula above).  Voxels which lie
outside this region will not be considered for filtering even if they lie
within the "-window" box.  Values between *1.5* and *3* are common.
The goal is to create a short-range filter which is as radially symmetric
as possible.
*This argument is optional and overrides the -window option.
If you specify either "-cutoff-xy" or "-cutoff-z", you must specify both.*

'''
   -cutoff-z  distanceZ
'''
specifies the width of the filter used in the Z direction (in units of "u" from the formula above, not voxels or *nm*).  Values between *1.5* and *3* are common.
*This argument is optional and overrides the -window option.
If you specify either "-cutoff-xy" or "-cutoff-z", you must specify both.*



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
