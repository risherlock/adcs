#include <stdio.h>
#include <math.h>

#include "attitude.h"
#include "time.h"
#include "iau06.h"

#define D2R 0.01745329251
#define R2D 57.2957795131
#define PI_2 1.57079632679489661923
#define PI 2 * 1.57079632679489661923

void print_dcm(const double dcm[3][3])
{
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      printf("  %.20f ", dcm[i][j]);
    }
    printf("\n");
  }
}

void frame_eci_to_ecef_dcm(const utc_t utc, double dcm[3][3])
{
  double W[3][3], R[3][3], Q[3][3];

  const double jd_since_j2000 = 51544.5;
  const double jd = time_julian_date(utc) - 2400000.5f;
  const double jd_elapsed = jd - jd_since_j2000;
  const double jd_fraction = fmod((fmod(jd_elapsed, 1.0) + 1.0), 1.0);
  double t = jd_elapsed / 36525.0;

  double xp = 0.0, yp = 0.0;
  const double s_prime = -0.000047 * t * D2R / 3600.0;
  const double zyx[3] = {s_prime, xp, yp};
  euler_to_dcm(zyx, EULER_ZYX, W);

  const double era = fmod(2.0 * PI * (jd_fraction + 0.7790572732640 + 0.00273781191135448 * jd_elapsed), 2.0 * PI);
  dcm_z(era, R);

  double x, y, s;
  iau06_get_xys(t, &x, &y, &s);

  double E = atan2(y, x);
  double d = atan(sqrt((x * x + y * y) / (1.0 - x * x - y * y)));
  const double zyz[3] = {E, d, -E - s};
  euler_to_dcm(zyz, EULER_ZYZ, Q);

  double WR[3][3];
  dcm_prod(W, R, WR);
  dcm_prod(WR, Q, dcm);
}

// void frame_ecef_to_ned(const double v_ecef[3], const double llr[3], double v_ned[3])
// {
//   const double st = sin(llr[0]);
//   const double ct = cos(llr[0]);
//   const double sp = sin(llr[1]);
//   const double cp = cos(llr[1]);
//   const double t  =  cp * v_ned[0] + sp * v_ned[1];

//   v_ned[0] = -st * t + ct * v_ned[2];
//   v_ned[1] = -sp * v_ned[0] + cp * v_ned[1];
//   v_ned[2] = -ct * t - st * v_ned[2];
// }

// void frame_ned_to_ecef(const double v_ned[3], const double llr[3], double v_ecef[3])
// {
//   const double st = sin(llr[0]);
//   const double ct = cos(llr[0]);
//   const double sp = sin(llr[1]);
//   const double cp = cos(llr[1]);
//   const double t  = - ct * v_ned[2] - st * v_ned[0];

//   v_ecef[0] = cp * t - sp * v_ned[1];
//   v_ecef[1] = sp * t + cp * v_ned[1];
//   v_ecef[2] = - st * v_ned[2] - ct * v_ned[0];
// }
