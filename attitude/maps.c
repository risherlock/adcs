#include <math.h>
#include "maps.h"

void quat_to_dcm(const double q[4], double r[3][3])
{
  r[0][0] = q[0] * q[0] + q[1] * q[1] - q[2] * q[2] - q[3] * q[3];
  r[0][1] = 2.0 * (q[1] * q[2] + q[0] * q[3]);
  r[0][2] = 2.0 * (q[1] * q[3] - q[0] * q[2]);
  r[1][0] = 2.0 * (q[1] * q[2] - q[0] * q[3]);
  r[1][1] = q[0] * q[0] - q[1] * q[1] + q[2] * q[2] - q[3] * q[3];
  r[1][2] = 2.0 * (q[2] * q[3] + q[0] * q[1]);
  r[2][0] = 2.0 * (q[1] * q[3] + q[0] * q[2]);
  r[2][1] = 2.0 * (q[2] * q[3] - q[0] * q[1]);
  r[2][2] = q[0] * q[0] - q[1] * q[1] - q[2] * q[2] + q[3] * q[3];
}

void dcm_to_quat(const double r[3][3], double q[4])
{
  double trace = r[0][0] + r[1][1] + r[2][2];
  double s;

  if (trace > 0)
  {
    s = 2.0 * sqrt(1.0 + trace);
    q[0] = 0.25 * s;
    q[1] = (r[2][1] - r[1][2]) / s;
    q[2] = (r[0][2] - r[2][0]) / s;
    q[3] = (r[1][0] - r[0][1]) / s;
  }
  else
  {
    if (r[0][0] > r[1][1] && r[0][0] > r[2][2])
    {
      s = 2.0 * sqrt(1.0 + r[0][0] - r[1][1] - r[2][2]);
      q[0] = (r[2][1] - r[1][2]) / s;
      q[1] = 0.25 * s;
      q[2] = (r[0][1] + r[1][0]) / s;
      q[3] = (r[0][2] + r[2][0]) / s;
    }
    else if (r[1][1] > r[2][2])
    {
      s = 2.0 * sqrt(1.0 + r[1][1] - r[0][0] - r[2][2]);
      q[0] = (r[0][2] - r[2][0]) / s;
      q[1] = (r[0][1] + r[1][0]) / s;
      q[2] = 0.25 * s;
      q[3] = (r[1][2] + r[2][1]) / s;
    }
    else
    {
      s = 2.0 * sqrt(1.0 + r[2][2] - r[0][0] - r[1][1]) * 2;
      q[0] = (r[1][0] - r[0][1]) / s;
      q[1] = (r[0][2] + r[2][0]) / s;
      q[2] = (r[1][2] + r[2][1]) / s;
      q[3] = 0.25 * s;
    }
  }
}

void quat_to_euler(const double q[4], double e[3], const euler_seq_t es)
{
  double q2_sq = q[2] * q[2];
  e[0] = atan2(2.0 * (q[0] * q[1] + q[2] * q[3]), 1.0 - 2.0 * (q[1] * q[1] + q2_sq));
  e[2] = atan2(2.0 * (q[0] * q[3] + q[1] * q[2]), 1.0 - 2.0 * (q2_sq + q[3] * q[3]));

  double sin_e1 = 2.0 * (q[0] * q[2] - q[1] * q[3]);

  if (fabs(sin_e1) >= 1.0)
  {
    e[1] = copysign(M_PI_2, sin_e1);
  }
  else
  {
    e[1] = asin(sin_e1);
  }
}

void euler_to_quat(const double e[3], const euler_seq_t es, double q[3])
{
  const double c[3] = {cos(0.5 * e[0]), cos(0.5 * e[1]), cos(0.5 * e[2])};
  const double s[3] = {sin(0.5 * e[0]), sin(0.5 * e[1]), sin(0.5 * e[2])};

  q[0] = c[0] * c[1] * c[2] + s[0] * s[1] * s[2];
  q[1] = c[0] * c[1] * s[2] - s[0] * s[1] * c[2];
  q[2] = s[0] * c[1] * s[2] + c[0] * s[1] * c[2];
  q[3] = s[0] * c[1] * c[2] - c[0] * s[1] * s[2];
}

void euler_to_dcm(const double e[3], const euler_seq_t es, double r[3][3])
{
  const double c[3] = {cos(e[0]), cos(e[1]), cos(e[2])};
  const double s[3] = {sin(e[0]), sin(e[1]), sin(e[2])};

  r[0][0] = c[0] * c[1];
  r[0][1] = s[0] * c[1];
  r[0][2] = -s[1];
  r[1][0] = c[0] * s[1] * s[2] - s[0] * c[2];
  r[1][1] = s[0] * s[1] * s[2] + c[0] * c[2];
  r[1][2] = c[1] * s[2];
  r[2][0] = c[0] * s[1] * c[2] + s[0] * s[2];
  r[2][1] = s[0] * s[1] * c[2] - c[0] * s[2];
  r[2][2] = c[1] * c[2];
}

void dcm_to_euler(const double r[3][3], double e[3], const euler_seq_t es)
{
  e[0] = atan2(r[0][1], r[0][0]);
  e[1] = asin(-r[0][2]);
  e[2] = atan2(r[1][2], r[2][2]);
}
