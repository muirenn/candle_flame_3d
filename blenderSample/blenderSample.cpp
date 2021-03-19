#include "FullAnimation.h"

int main()
{
	const int T = N_FRAMES;
    const int n = M_TO_VXL(N_SIDE);
    cout << T << "  " << n << endl;
	const string filename = "sample20212.bvox";

	time_t start_time;
	struct tm * timeinfo;
	time(&start_time);
	timeinfo = localtime(&start_time);

	char log_filename[LOG_NAME_LENGTH];
	strftime(log_filename,LOG_NAME_LENGTH,"%d%b%Y_%H%M%S.log",timeinfo); // DayMonthYear_HourMinuteSecond.log
	cout << log_filename << endl;

	FullAnimation* animation = new FullAnimation(T, n, filename, log_filename);
	animation->init();
	animation->run();
	animation->write_animation_to_file();

	delete animation;
	return 0;
}
 