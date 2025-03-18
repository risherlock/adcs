#include <stdio.h>
#include <math.h>

#include "attitude.h"
#include "time.h"
#include "iau06.h"

#define D2R 0.01745329251
#define R2D 57.2957795131
#define PI_2 1.57079632679489661923
#define PI 2 * 1.57079632679489661923
#define FP 1 / 298.257223563
#define RE 6378137
#define E2 FP * (2 - FP)

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
  const double euler_zyx[3] = {s_prime, xp, yp};
  euler_to_dcm(euler_zyx, EULER_ZYX, W);

  const double era = fmod(2.0 * PI * (jd_fraction + 0.7790572732640 + 0.00273781191135448 * jd_elapsed), 2.0 * PI);
  dcm_z(era, R);

  double x, y, s;
  iau06_get_xys(t, &x, &y, &s);

  double E = atan2(y, x);
  double d = atan(sqrt((x * x + y * y) / (1.0 - x * x - y * y)));
  const double euler_zyz[3] = {E, d, -E - s};
  euler_to_dcm(euler_zyz, EULER_ZYZ, Q);

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
double hypot(double x, double y) {
  return sqrt(x * x + y * y);
}

// Function to compute sign
double sign(double x) {
  return (x > 0) ? 1.0 : ((x < 0) ? -1.0 : 0.0);
}

// Function to convert ECEF to geodetic coordinates
void frame_ecef_to_lla(const double ecef[3],double lla[3])
{
  const double x = ecef[0];
  const double y = ecef[1];
  const double z = ecef[2];

  const double rho = hypot(x,y);

  const double a = RE;         // Semimajor axis (meters)
  const double f = FP; // Flattening
  const double b = (1 - f) * a;       // Semiminor axis
  const double e2 = E2;      // Square of (first) eccentricity
  const double ae2 = a * e2;          // a * e^2
  const double bep2 = b * e2 / (1 - e2);

  // Calculate longitude (lambda)
  lla[1] = atan2(y, x);

  // Starting value for parametric latitude (beta)
  double r = hypot(rho, z);
  double u = a * rho;
  double v = b * z * (1 + bep2 / r);

  double cosbeta = sign(u) / hypot(1, v / u);
  double sinbeta = sign(v) / hypot(1, u / v);

  // Fixed-point iteration with Bowring's formula
  int count = 0;
  double cosprev, sinprev;

  do {
    cosprev = cosbeta;
    sinprev = sinbeta;

    // Update u and v
    u = rho - ae2 * pow(cosbeta, 3);
    v = z + bep2 * pow(sinbeta, 3);

    // Update cosbeta and sinbeta
    double au = a * u;
    double bv = b * v;
    cosbeta = sign(au) / hypot(1, bv / au);
    sinbeta = sign(bv) / hypot(1, au / bv);

    // Check for convergence
    count++;
  } while((hypot(cosbeta - cosprev, sinbeta - sinprev) > 1e-10) && (count < 5));

  // Calculate geodetic latitude
  lla[0] = atan2(z + bep2 * pow(sinbeta, 3), rho - ae2 * pow(cosbeta, 3));

  // Calculate altitude 
  double N = a / sqrt(1 - e2 * pow(sin(lla[0]), 2));
  lla[2] = rho * cos(lla[0]) + (z + e2 * N * sin(lla[0])) * sin(lla[0]) - N;

  lla[0] *= R2D;
  lla[1] *= R2D;
}
