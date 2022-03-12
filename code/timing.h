#ifndef TIMING_H
#define TIMING_H

#include <time.h>

clock_t ticc_314159265358979323846;
time_t tict_314159265358979323846;

void tic()
  {
    ticc_314159265358979323846=clock();
  }

void toc()
  {
    float ts= (double)(clock() - ticc_314159265358979323846) / CLOCKS_PER_SEC;
    printf("Elapsed time: %f s\n",ts);
  }

  double toq()
  {
    return (double)(clock() - ticc_314159265358979323846) / CLOCKS_PER_SEC;
  }

#endif