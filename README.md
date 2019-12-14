# candle_flame_3d

To visualize 3D data, it is necessary to define a finite volume
within the Blender scene; this is done by creating a cube object
that will be the domain for the visualization. Once the voxel
data file is assigned to this domain, the 3D information stored
in the file will be rendered in the Blender scene within the
boundaries of the cube. 

The voxel file represents a datacube subdivided in the three
cartesian axes (x, y, z). For each discrete cartesian coordinate (i,
j, k ) inside the volume, there is one voxel value V(i, j, k).
Additionally, it is possible to include multiple snapshots
(frames) in the same file. This allows showing the evolution of
the datacube in time. Therefore, the voxel data is defined by
• The resolution numbers (nx, ny, nz), which indicate the
  number of subdivisions of the domain in each axis;
• The frame number nf, which indicates the total number of
  frames contained in the file; and
• (nx × ny × nz) × nf values, which describe the data that will be 
  visualized inside the domain.

The Blender Voxel format consists in a binary file that has
4 × 32bit integers in the header, corresponding to (nx, ny, nz,
nf), followed by the (nx × ny × nz × nf) × 32bit floating point
values that represent the actual data. For this format, the values
must be normalized to a [0, 1] interval. The order to write the
full datacube after the header is
1. Frame by frame, where each frame has nz layers.
2. Layer by layer, where each layer has ny lines.
3. Line by line, where each line has nx values.
4. Value by value, until completing the nx values.


https://iopscience.iop.org/article/10.1088/1538-3873/129/975/058010/pdf


Normalization: zi = (xi - min(x))/(max(x)-min(x));
