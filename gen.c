#include <stdlib.h>
#include <time.h>
#include <stdio.h>

int main() {
	int c = 10000000;
	int i;

	//initstate(time(NULL), "chnyangzhe-programming@yahoo.com.cn", 32);
	srandom(time(NULL));

	for (i = 0; i < c; ++i)
		putchar('a' + ((random() >> 2) & 15));
	puts("");
}
