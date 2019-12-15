#include "stdafx.h"
#include "FullAnimation.h"

int main()
{
	const int T = 30;
	const int n = 50;
	const string filename = "holyshit4.bvox";

	FullAnimation* animation = new FullAnimation(T, n, filename);
	animation->init();
	animation->run();
	animation->write_animation_to_file();

	delete animation;
	return 0;
}
