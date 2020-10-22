#include "FullAnimation.h"

int main()
{
	const int T = N_FRAMES;
    const int n = M_TO_VXL(N_SIDE);
    cout << T << "  " << n << endl;
	const string filename = "sample9.bvox";

	time_t start_time;
	struct tm * timeinfo;
	time(&start_time);
	timeinfo = localtime(&start_time);

	char log_filename[20];
	strftime(log_filename,20,"%d%b%C_%H%M%S.log",timeinfo);
	cout << log_filename << endl;

	FullAnimation* animation = new FullAnimation(T, n, filename, log_filename);
	animation->init();
	animation->run();
	animation->write_animation_to_file();

	delete animation;
	return 0;
}
 