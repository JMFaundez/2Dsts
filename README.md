# Average statistics 2D
Used after running airfoil simulation and collect 2d stats. Main goal is to get the data for turbulence intesity and friction coefficient.

## MATLAB
### Jsts.m
It selects the sts files to average and call `interpolate_data.m` for each file to interpolate the data over the right mesh.

### interpolate_data.m
Read the right 2d mesh and interpolate the corresponding sts file over the mesh

### data_to_nek.m 
Read the interpolated sts file and write it as Nek field.

### matlab_script
Folder that contains some useful code related to nek

### plot_sts.m
Plot stuff after havinf the data interpolated and the gradients computed by nek.

### TI_data.m
Read the interpolated sts and write the useful information used to compute the turbulence intensity

## NEK
The nek files are used to compute the gradient of u' and w' (instead of using finite differences in matlab). It needs the field with the cosine and sine of the transformation to the natural coordinates.


