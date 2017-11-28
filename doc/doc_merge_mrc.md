
merge_mrc
===========
**merge_mrc** is a program for combining two volumetric images (i.e. tomograms, both of identical size) into one image/tomogram, using a combination of addition, multiplication, and thresholding operations.  These features can be used perform binary operations between two images which are similar to "**and**" and "**or**" operations.  ("**not**" operations are also possible.  See below.) As with the "filter_gauss" and "filter_dog" programs, you can also use the "-mask" argument to restrict the operation to certain voxels from the image.
*(A detailed description of the threshold and rescaling functions used in this step are provided at the end of the the "doc_filter_gauss.md" file.)*
Detailed documentation for this program is located in in the "*doc/*" subdirectory.

### Clipping the output range
*Note:*  After adding or multiplying the voxel brightnesses, the brightness of each resulting voxel is automatically clipped between 0 and 1, *before* saving the result to a file.
(In other words, the resulting brightnesses below 0 are replaced with 0, and the resulting brightnesses above 1 are replaced with 1.
*Perhaps I will provide a more general way to rescale and clip the output brightnesses later...*)


### merge_mrc examples:
To **add** the brightness values between to tomograms and save the result in a new file ("out_file.mrc"), use:
```
   merge_mrc file1.mrc  +  file2.mrc  out_file.mrc
```
To **multiply** the brightness values use:
```
   merge_mrc file1.mrc "*" file2.mrc  out_file.mrc
```
The quotes around the ** * ** character are

The following examples apply 1, 2 or 4 thresholds to the input tomograms before performing the **+** or **\*** operation.  You specify the thresholds you want to use with commas placed after the file name (no spaces).  The following command will replace all of the brightness of all the voxels whos brightess is below 0.5 with 0, and all of brightnesses above 0.5 with 1, and *then* multiply them together:
```
   merge_mrc file1.mrc,0.5 "*" file2.mrc,0.5 out_file.mrc
```
(Note that 8-bit and 16-bit integer brightnesses are replaced with floating point numbers in the range from 0 to 1 beforehand, and the resulting tomogram is saved in 32bit float format.)


### Applying threshold filters to the input images:

It is often convenient to rescale or clip the brightnesses of the voxels from either image (tomogram) *before* adding or multiplying them together.  This way, you can insure that most of the voxels in the image are either ***0*** or ***1*** beforehand.  (This way, multiplying and adding the resulting voxel brightnesses is equivalent to performing an "and" and "or" gate operations on these 0,1 values.)  The resulting image (tomgram) created by *merge_mrc* can be useful for segmentation (or masking).

The following command will replace all of the brightness of all the voxels whose brightess falls below *0.48* with ***0***, and all of brightnesses above *0.49* with ***1*** (linearly scaling any voxels with brightnesses between *0.48* and *0.49* to fill the range from *0* to *1*).  *Then* it will add them together
```
   merge_mrc file1.mrc,0.48,0.49 + file2.mrc,0.48,0.49 out_file.mrc
```
*(Again, a detailed description of the threshold and rescaling functions used in this step are provided at the end of the the "doc_filter_gauss.md" file.)*

When 4 numbers follow an input file name, the brightness of all the voxels from that tomogram in that file will be run through a **double-threshold** filter.  For example to replace voxels from both tomograms whose brightess falls below *0.48* with ***0***, and replace all of brightnesses between *0.49* and *0.51* with ***1***, and all of the brightness values above *0.52* with ***0*** again (linearly scaling any voxels with brightnesses between *0.48* and *0.49*, or between *0.51* and *0.52*), use this command:
```
   merge_mrc file1.mrc,0.48,0.49,0.51,0.52 + file2.mrc,0.48,0.49,0.51,0.52 out_file.mrc
```
resulting in:
```
 output
 brightness
 (a.k.a. "density")
  /|\
 1 |                 ________________                
   |             _,-'                `-._
   |         _,-'                        `-._
 0 |______,-'                                `-._______\ input
        0.48       0.49            0.51       0.52     / brightness
```

*Note:* If the order of thresholds is reversed, the inverse images is generated.  For example, the following command will replace all of the brightness of all the voxels whose brightess falls below 0.4 with 1, and all of brightnesses above 0.6 with 0 (linearly scaling any voxels with brightnesses between 0.4 and 0.6 to fill the range from 1 to 0).  (If the original tomogram had brightnesses of either 1 or 0, this can be used to perform something analogous to a "not" operation.)  *Then* it will add the two tomograms together
```
   merge_mrc file1.mrc,0.6,0.4 + file2.mrc,0.6,0.4 out_file.mrc
```
Similarly, when 4 thresholds are specified, you can invert the output image by altering their order:
```
   merge_mrc file1.mrc,0.5,0.6,0.3,0.4 + file2.mrc,0.5,0.6,0.3,0.4 out_file.mrc
```
...which results in:
```
  /|\                                                   
 1 |_____                                       _______\ input
   |     `-._                               _.-'       / brightness
   |         `-._                       _,-'             ("density")
 0 |             `-._________________,-'
        0.3       0.4               0.5       0.6
```

### "not" operations
*In addition to "**and**" and "**or**" operations, you can also perform "not" operations.  You can perform the "**not**" operation on each **input** by reversing the order of the thresholds following a file name.  The example above ("file1.mrc,0.6,0.4) demonstrates how to do that.*

*"**not**" operations can also be performed on the **output** image by saving it to a file, and later using the "filter_gauss" program with a small sigma parameter of approximately ~1e-06 (to avoid blurring the image), and then using the "-thresh2" argument with the first threshold number greater than the second.*
