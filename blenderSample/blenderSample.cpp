#include "FullAnimation.h"

int main()
{
	const int T = N_FRAMES;
	const int n = M_TO_VXL(N_SIDE);
	cout << T << " frames;  " << n << " voxels per side" << endl;

	FS::path *log_filename = new FS::path(), *bvox_filename = new FS::path();

	get_filenames(log_filename, bvox_filename);
	cout << log_filename->string() << endl;

	FullAnimation *animation = new FullAnimation(T, n, *bvox_filename, *log_filename);
	animation->init();
	animation->run();
	animation->write_animation_to_file();
	// animation->write_fake_animation_to_file();

	delete animation;
	return 0;
}
