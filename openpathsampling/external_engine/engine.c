#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(int argc, char ** argv)
{
    if (argc < 2) {
        printf("Requires two arguments: delay time (ms) and filename");
        exit(1);
    }
    // argv[0] is program name
    int milliseconds = atoi(argv[1]);
    FILE * f = (argc == 2) ? stdout : fopen(argv[2], "w");

    struct timespec sleep_time, foo;

    sleep_time.tv_sec = milliseconds / 1000;
    sleep_time.tv_nsec = (milliseconds % 1000) * 1000000;

    int max_steps = 100000;
    int i; for (i=0; i<max_steps; i++) {
        fprintf(f, "0.0\n");
        nanosleep(&sleep_time, &foo);
    }
}
