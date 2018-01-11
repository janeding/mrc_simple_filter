filter_mrc
===========

**filter_mrc** applies a filter to a tomogram in the X,Y,Z directions
and saves the result as a new .mrc/.rec file.
This program can be used to rescale or invert a 3-D image,
remove high or low frequencies (smoothing, edge detection, band pass filter).
perform 3-D blob detection.
Currently, the program supports the following list of filters:
(generalized) Gaussians,
(generalized) Difference-of-Gaussians,
Laplacian-of-Gaussians, and others.
Both isotropic and anisotropic filters are supported.


This program supports masks and thresholding using
the "*-mask*" and "*-thresh*" (and "*-thresh2*" and "*-thresh4*") arguments.

*In this program all of the parameters provided by the user are assumed to be in **physical units**.  (Not **voxels**.  This effects all of the gaussian "sigma" parameters, specifically all of the "-width-xy" and "-width-z" and "-window-ratio" arguments.  To express these distances in *voxels* instead, set the voxel_width to 1.0 using "-w 1.0".  See below)*


## Usage Example:

```
filter_mrc \
   -in Bdellovibrio_1K.rec \
   -out Bdellovibrio_1K_ribosome_peaks.mrc \
   -w 1.92 \
   -dog 12.0 12.0 24.0 15.0 15.0 30.0 \
   -exponents 2 2 \
   -cutoff 0.01 \
   -mask Bdellovibrio_1K_mask_water=0_periplasm=1_cytoplasm=2.mrc \
   -mask-select 0 -mask-out 1.0
```

## Arguments:

### Input and Output files

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

### Voxel Width

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

The "**-gauss**" filter uses a (low-pass)
[Gaussian blur](https://en.wikipedia.org/wiki/Gaussian_blur)
filter.
(The "**-gauss-aniso**" filter allows you to apply different amounts of blurring to the X,Y, and Z directions.)

The "**-dog**"" filter uses a (band-pass)
[Difference-Of-Gaussians](https://en.wikipedia.org/wiki/Difference_of_Gaussians)
filter
which can be used for blob detection as well as low-frequency removal.
(The "**-dog-aniso**" filter allows you to modify the properties of the filter in the X,Y,Z directions.)

The "**-log**" filter can be used for [blob detection](https://en.wikipedia.org/wiki/Blob_detection#The_Laplacian_of_Gaussian).

In all of these filters, [*Generalized Gaussians*](https://en.wikipedia.org/wiki/Generalized_normal_distribution#Version_1)
can be substituted for ordinary Gaussians by using the "**-exponent m**" or "**-exponents m n**" command line arguments.  (Default value 2.  See below.)  Increasing these exponents sharpens the edges of the filter

### Template Matching

You can Difference-Of-Generalized-Gaussians filter to create a filter 
whose value is a positive constant within a spherical volume, 
and is negative outside this volume within a larger radius. 
("**-dog**" and "**-exponents**") filter can be used for the template matching of compact (approximately spherical) objects.  Such filters can be used to identify simple objects such as ribosomes or nucleosomes according to their approximate size.

However generalized Gaussians are orders of magnitude slower than ordinary Gaussian filters.  For this reason a hybrid filter, "**-dogxy**", is available which uses a Difference-Of-Generalized-Gaussians filter in the XY direction combined with an ordinary Gaussian filter in the Z direction.  The computational cost of "**-dogxy**" lies in between "**-dog**" and "**-dog**".  (Features in electron tomography are typically blurred more in the Z direction due to the effect of the missing wedge.  So it may be pointless and impossible to use the computationally more expensive "**-dog**" filter in an effort to find the precise boundaries of objects in the Z direction.  In these cases, the "**-dogxy**" filter may work just as well and is significantly faster.)


### -gauss
The **-gauss** and **-gauss-aniso** arguments must be followed by one or more numbers specifying the width of the Gaussian filter to apply to your image:
```
   -gauss-aniso  s_x  s_y  s_z
```
```
   -gauss  s      (this means s_x = s_y = s_z = s)
```
When the "-gauss" or "-gauss-aniso" filter is selected, the
original image is convolved with the following function:
```
   h(x,y,z) = A*exp(-0.5 * r^2)
    where r = sqrt((x/s_x)^2 + (y/s_y)^2 + (z/s_z)^2))
```
The width of the Gaussian (ie, the s_x, s_y, s_z arguments) should be specified in units of physical distance, not in voxels.
(The A coefficient will determined automatically by normalization, ie., so that the discrete sum of h(x,y,z) over x,y,z is 1.)
The width of the filter-window is chosen automatically according the s, s_x, s_y, s_z parameters selected by the user.  However this can be customized using the "-cutoff" and "-window-ratio" arguments if necessary (see below).

*2018-1-04: The feature described below has not yet been implemented, but is trivial to add:*

*Note:* If the optional "**-exponent m**" argument is supplied, then
the original image is instead convolved with
[the following function](https://en.wikipedia.org/wiki/Generalized_normal_distribution#Version_1):
```
h(x,y,z) = A*exp(-r^m)
 where r = sqrt((x/s_x)^2 + (y/s_y)^2 + (z/s_z)^2))
```
This crude generalization of the "gaussian" function gives us an ad-hoc way
to alter shape of the filter.
(As you increase the **n** exponent parameter,
the shape of the filter becomes sharper
and begins to resemble a step-function.
For your convenience, you can plot these functions for various parameters using
the "draw_filter1D.py -ggauss A a m" python script which is located in the same directory.
*Warning: This will slow down the filter considerably* because in this case the filter is no longer a [separable](https://en.wikipedia.org/wiki/Separable_filter) function of x,y,z.)


### -log
The **-log** and **-log-aniso** argument must be numbers specifying the approximate radius (width/2) of features that you wish to emphasize in the image:

```
  -log-aniso  s_x  s_y  s_z
```
```
  -log  s      (this means s_x = s_y = s_z = s)
```
When the "**-log**" or "**-log-aniso**" filter is selected,
a ["laplacian-of-gaussians"](https://en.wikipedia.org/wiki/Blob_detection#The_Laplacian_of_Gaussian)
 filter will be used.
The original image is convolved with the following function
(which is an approximation to the "laplacian-of-gaussian" filter):
```
   h(x,y,z) = A*exp(-0.5 * r_a^2) - B*exp(-0.5 * r_b^2)
  where r_a = sqrt((x/a_x)^2 + (y/a_y)^2 + (z/a_z)^2))
    and r_b = sqrt((x/b_x)^2 + (y/b_y)^2 + (z/b_z)^2))
      a_x=s_x*(1-0.5*delta), a_y=s_y*(1-0.5*delta), a_z=s_z*(1-0.5*delta)
      b_x=s_x*(1+0.5*delta), b_y=s_y*(1+0.5*delta), b_z=s_z*(1+0.5*delta)
```
The A and B parameters are determined automatically by normalization.
The "delta" parameter is 0.01 by default.
(This can be overridden using the "-log-delta delta" argument.
A smaller "delta" value may improve the approximation.)
The width of the Gaussian (the s_x, s_y, s_z arguments) should be specified in units of physical distance, not in voxels.
The A and B coefficients will be automatically chosen so that the discrete sum of h(x,y,z) over x,y,z is 0, and the peak height is 1 (A-B=1).
The width of the filter-window is chosen automatically according the shape of the filter selected by the user.  However this can be customized using the "-cutoff" and "-window-ratio" arguments if necessary (see below).


### -dog-aniso
The **-dog** and **-dog-aniso** arguments must be followed by numbers indicating the width of the Gaussians you wish to convolve with your image.  (Equivalently, the numbers indicate the band of frequencies you wish to keep from your original image.)

```
  -dog-aniso  a_x  a_y  a_z  b_x  b_y  b_z
```
```
  -dog  a  b      (this means a_x=a_y=a_z=a and b_x=b_y=b_z=b)
```
When the "**-dog**" or "**-dog-aniso**" filter is selected, a
[Difference-of-Gaussians](https://en.wikipedia.org/wiki/Difference_of_Gaussians)
filter will be used.
The original image will be convolved with the following function:
```
   h(x,y,z) = A*exp(-0.5 * r_a^2) - B*exp(-0.5 * r_b^2)
  where r_a = sqrt((x/a_x)^2 + (y/a_y)^2 + (z/a_z)^2))
    and r_b = sqrt((x/b_x)^2 + (y/b_y)^2 + (z/b_z)^2))
```
The width of the Gaussian (the a_x, a_y, a_z, b_x, b_y, b_z arguments) should be specified in units of physical distance, not in voxels.
The A and B coefficients will be automatically chosen so that the discrete sum of h(x,y,z) over x,y,z is 0, and the peak height is 1 (A-B=1).
The width of the filter-window is chosen automatically according the shape of the filter selected by the user.  However this can be customized using the "**-cutoff**" and "**-window-ratio**" arguments if necessary (see below).

*Note:* If the optional "**-exponents m n**" argument is supplied, then
the original image is instead convolved with the following function:
```
h(x,y,z) = A*exp(-r_a^m) - B*exp(-r_b^n)
 where r_a = sqrt((x/a_x)^2 + (y/a_y)^2 + (z/a_z)^2))
   and r_b = sqrt((x/b_x)^2 + (y/b_y)^2 + (z/b_z)^2))
```
This crude generalization of the "gaussian" function gives us an ad-hoc way
to alter shape of the filter.
(As you increase the **m** and **n** exponent parameters,
the shape of the filter becomes sharper
and begins to resemble step-functions.
For your convenience, you can plot these functions for various parameters using
the "draw_filter1D.py -dogg A B a b m n" python script which is located in the same directory.
*Warning: This will slow down the filter considerably* because in this case
the filter is no longer a
[separable](https://en.wikipedia.org/wiki/Separable_filter) function of x,y,z.)


### -dogxy
The **-dogxy** argument must be followed by 3 numbers:
```
  -dogxy  a  b  c
```
When the "**-dogxy**" filter is selected,
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
 using the "**-exponents m n**" argument.

 Along the Z direction, the filter used is a simple Gaussian:
```
   h_z(r) = C * exp(-0.5*(z/c)^2)
```
(The "C" constant is determined by normalization.)



## Thresholding and Rescaling (Intensity Map Filters):

Thresholding filters provide a way to rescale (map) the 
intensity of all of the voxels in the image to new intensities.
For example, you could multiply these intensities by a constant, 
and/or clip them within a specified range.
Thresholding operations can be used to insure that the brightnesses 
in the image lie between 0 and 1.
They can also be used to select only voxels whose intensities lie within 
a narrow range of interest.

(These operations are quick to calculate because each
 voxel's intensity depends only on its previous intensity,
 *not* the intensities of any of the voxels nearby.)

*(Again, voxels from 8-bit or 16-bit integer images are converted to floating point numbers between 0 and 1 beforehand, by dividing them by 255 and 65535, respectively.  Keep this in mind when using software such as IMOD which typically reports pixel brightnesses between 0 and 255.)*

*Note:* All threshold operations are performed *after* any other filtering operations have been applied (such as Gaussian blurs, -dog filters, or -log filters).


### -invert

This filter replaces bright voxels with dark ones, and visa versa, 
while keeping the average voxel brightness the same.
(To calculate the new voxel intensity, 
 the *difference* between the original voxel intensity and the average intensity
 is *subtracted* from the average intensity.)
Inverting an image can be useful to prepare files which you plan 
to read with other software (such as ChimeraX and IMOD) 
which interpret or display voxel intensities differently.


### -thresh  thresh01

If the "-thresh" argument is passed, it must be followed by a number ("thresh01").  Voxels with intensities *below* this number will be replaced 0, and voxels *above* this number will be replaced with 1.  For example:
```
-thresh 0.5
```
results in:
```
 output
 intensity
  /|\                      _____________________________\
 1 |                      |                             /
   |                      |                 
   |                      |            
 0 |______________________|  __________________________\ input
                          ^                            / intensity
                         0.5
                     (threshold)
                                                          
```

### -thresh2  thresh_a  thresh_b

If the "-thresh2" argument is passed, then it must be followed by 2 numbers.
For example:
```
-thresh2 0.48 0.52
```
In this case, the resulting voxel intensities in the range 
from *thresh_a* to *thresh_b* will be scaled between 0 and 1 
and clipped above and below.
Graphically, the threshold filter resembles a step function
with a smooth ramp between *thresh_a* and *thresh_b*:

```
 output
 intensity
  /|\                              _________________
 1 |                           _.-'                 
   |                       _,-'                 
   |                   _,-'            
 0 |________________,-'                     ________\ input
                    ^              ^                / intensity
                   0.48          0.52
                (thresh_a)     (thresh_b)
```

***Alternatively***, if the threshold parameters are reversed
(if ***thresh_b < thresh_a***), as in this example:
```
-thresh2 0.52 0.48
```
...then the output intensity will be inverted
(i.e., bright voxels become dark, and dark voxels become bright):

```
 output
 intensity
  /|\
   |________________
 1 |                `-._
   |                    `-._
   |                        '-._            
 0 |                            `-.___________________\ input
                    ^              ^                  / intensity
                   0.48          0.52
                (thresh_b)     (thresh_a)
```


### -thresh4  thresh01_a  thresh01_b  thresh10_a  thresh10_b

Sometimes it is useful to select a ***narrow range of voxel intensities***.
*filter_mrc* allows you to do this using the **"-thresh4"** argument.
The **"-thresh4"** must be followed by 4 numbers in order, for example:
``` 
 -thresh4 0.3 0.4 0.5 0.6
```
The resulting image voxels will be scaled
between 0 and 1 according to the following function:

```
 output
 intensity
  /|\
 1 |                 ________________                
   |             _,-'                `-._
   |         _,-'                        `-._
 0 |______,-'                                `-._______\ input
         ^         ^                 ^         ^       / intensity
        0.3       0.4               0.5       0.6
   thresh_01_a  thresh_01_b    thresh_10_a  thresh_10_b
```
This assumes the numbers were listed in *increasing* order.
***Alternatively***, if the parameters are listed in *decreasing* order,
for example:
``` 
 -thresh4 0.6 0.5 0.4 0.3
```
...then the output is inverted:
```
 output
 intensity
  /|\                                                   
 1 |_____                                       _________
   |     `-._                               _.-'       
   |         `-._                       _,-'             
 0 |             `-._________________,-'          _______\ input
         ^         ^                 ^         ^         / intensity
        0.3       0.4               0.5       0.6
   thresh_10_b  thresh_10_a    thresh_10_b  thresh_10_a
```



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




## Additional Arguments:

### Filter Size
```
   -cutoff threshold
```
This specifies the number of voxels in the filter which will be convolved
with the image.
The filter window width is extended in the X, Y, and Z directions until
the filter intensity decays to less than **threshold** (a fraction),
compared to its peak value.
(For Difference-of-Gaussian filters, both Gaussians must decay to this value.)
If unspecified, this parameter is set to 0.02.
(This overrides the "-window-ratio" argument.)

```
   -window-ratio Wx Wy Wz
```
If you prefer to specify the size of the filter window manually, you can
use the "-window-ratio" argument.
This argument specifies the size of the 3-D filter box in x,y, directions
(in units of the **a** and **b** distance parameters, or **s**(sigma)
 parameters which appear in the formulas above).
This overrides the "-cutoff threshold" argument.

*(Note: Incidentally, the product of these 3 numbers, Wx*Wy*Wz, is proportional to the running time of the filter.
However when using the "-gauss" or "-dog" filters with default exponent settings,
the running time is proportional to the sum of these numbers, Wx+Wy+Wz.
Keep this in mind when specifying filter window widths.)*


```
   -ang-to-nm
```
When reading MRC files, this argument will
multiply the voxel width
specified in the MRC file by 0.1.
This allows you to specify your distance parameters in units of nm
instead of Angstroms.

This argument has no effect if the "-w" argument is used.



