#include <stdio.h>
#include "frame.h"
#include "time.h"

int main()
{
  double dcm[3][3];
  utc_t t = {.year = 2000, .month = 1, .day = 12, .hour = 4, .minute = 52, .second = 12};
  frame_eci_to_ecef_dcm(t, dcm);

  return 0;
}
