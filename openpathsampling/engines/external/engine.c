#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/*
 * engine.c
 *
 * Trivial 1D "engine" to be used for tests of an external engine.
 *
 * Allows us to control pace of creating new frames, as well as initial
 * frame condition.
 */

int main(int argc, char ** argv)
{
    if (argc < 2) {
        printf("Requires two arguments: delay time (ms) and filename\n");
        exit(1);
    }
    // argv[0] is program name
    int milliseconds = atoi(argv[1]);
    FILE * f = (argc == 2) ? stdout : fopen(argv[2], "w");
    double initial_position = 0.0;
    double velocity = 1.0;
    if (argc == 4) {
        // this means we have an input file given last
        FILE * in_f = fopen(argv[3], "r");
        fscanf(in_f, "%lf %lf", &initial_position, &velocity);
    }

    struct timespec sleep_time, foo;

    sleep_time.tv_sec = milliseconds / 1000;
    sleep_time.tv_nsec = (milliseconds % 1000) * 1000000;

    int max_steps = 1000000;
    double position = initial_position;
    int i; for (i=0; i<max_steps; i++) {
        position += velocity;
        fprintf(f, "%lf %lf\n", position, velocity); fflush(f); 
        nanosleep(&sleep_time, &foo);
    }
}
