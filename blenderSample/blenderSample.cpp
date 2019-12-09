#include "stdafx.h"

/*
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
*/

// #include "VolumeFrame.h"
#include "VolumeInstance.h"

constexpr auto T = 20;;

int main()
{
	// RenderedVolume sample_volume;

	// sample_volume.fill_in_volume_sample1();

	// RenderedVolume test("testingstam.bvox",100, 1);
	/* VolumeFrame* frame = new VolumeFrame();

	//float *density, *density_prev, *velocity_x_prev, *velocity_y_prev, *velocity_z_prev,
	//*velocity_x, *velocity_y, *velocity_z;
	//float visc = 0.00001, dt = 0.1, diff = 0;


	frame->arbitrary_setup();
	frame->run();
	
	delete frame;
	*/

	VolumeInstance* volume_instance = new VolumeInstance();

	volume_instance->init();

	ofstream blender_file;
	blender_file.open(volume_instance->filename_, ios::out | ios::binary | ios::app);
	blender_file.write(reinterpret_cast<const char*>(&T), sizeof(T));
	blender_file.close();

	volume_instance->write_to_file(volume_instance->density_);

	for (int i = 1; i < T; i++)
	{
		volume_instance->step();
	}

	delete volume_instance;
	return 0;
}
