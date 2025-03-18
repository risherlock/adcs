// Use case of IGRF-14 implementation.
// 2021-12-17

#include <stdio.h>
#include "igrf.h"

int main()
{
  // Miss Violet Smith disturbed Holmes 130 years ago on this day
  const igrf_time_t dt = {.year = 2025, .month = 04, .day = 23, 0, 0, 0};

  // Geodetic coordinates of 221B Baker Street on LEO
  const double latitude = 51.523788; // deg
  const double longitude = -0.158611; // deg
  const double altitude = 400.0; // km
  const double x[3] = {latitude, longitude, altitude};

  // Magnetic field in NED frame
  double b[3] = {0.0};
  bool status = igrf(dt, x, IGRF_GEODETIC, b);

  if (status)
  {
    printf("Inputs:\n");
    printf("  Time      : %d-%d-%d, %d:%d:%d\n", dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second);
    printf("  Latitude  : %f deg\n", latitude);
    printf("  Longitude : %f deg\n", longitude);
    printf("  Altitude  : %f km\n", altitude);
    printf("\nOutputs:\n");
    printf("  Bn          : %f nT\n", b[0]);
    printf("  Be          : %f nT\n", b[1]);
    printf("  Bd          : %f nT\n", b[2]);
    printf("  Magnitude   : %f nT\n", igrf_mag(b));
    printf("  Inclination : %f deg\n", igrf_inc(b) * R2D);
    printf("  Declination : %f deg\n", igrf_dec(b) * R2D);
  }
  else
  {
    printf("Date error!\n");
  }

  return 0;
}

/*
Expected output:

Inputs:
  Time      : 2025-4-23, 0:0:0
  Latitude  : 51.523788 deg
  Longitude : -0.158611 deg
  Altitude  : 400.000000 km

Outputs:
  Bn          : 16625.026351 nT
  Be          : 69.051911 nT
  Bd          : 37532.201921 nT
  Magnitude   : 41049.512182 nT
  Inclination : 66.108694 deg
  Declination : 0.237976 deg
*/
