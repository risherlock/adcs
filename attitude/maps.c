#include <math.h>
#include <stdbool.h>
#include <inttypes.h>

#include "maps.h"

#define M_PI_2 1.57079632679489661923
#define M_PI 2 * 1.57079632679489661923

// Choose an angle on infinite solutions
#define ANGLE_ON_SINGULARITY_RAD 0.0f

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

// Quaternion to intrinsic Euler angles as described in Ref.[2].
void quat_to_euler(const double q[4], double e[3], const euler_seq_t es)
{
  const double tolerance = 1e-7;

  // Parse Euler sequence
  const uint8_t i = (es / 100) % 10;
  const uint8_t j = (es / 10) % 10;
  uint8_t k = es % 10;

  // Tait-Bryan angles
  bool not_proper = true;

  // Proper Euler angles
  if (i == k)
  {
    k = 6 - i - j;
    not_proper = false;
  }

  // Is permutation even or odd?
  const int epsilon = -(i - j) * (j - k) * (k - i) / 2.0f;
  double a, b, c, d;

  if (not_proper)
  {
    a = q[0] - q[j];
    b = q[i] + q[k] * epsilon;
    c = q[j] + q[0];
    d = q[k] * epsilon - q[i];
  }
  else
  {
    a = q[0];
    b = q[i];
    c = q[j];
    d = q[k] * epsilon;
  }

  const double a_sq = a * a;
  const double b_sq = b * b;
  const double c_sq = c * c;
  const double d_sq = d * d;
  const double hyp_ab = a_sq + b_sq;
  const double hyp_cd = c_sq + d_sq;

  e[1] = acos(2.0 * ((hyp_ab) / (hyp_ab + hyp_cd)) - 1.0);
  const double theta_plus = atan2(b, a);
  const double theta_minus = atan2(d, c);

  // Check singularity
  if (fabs(e[1]) < tolerance)
  {
    e[0] = ANGLE_ON_SINGULARITY_RAD;
    e[2] = 2 * theta_plus - ANGLE_ON_SINGULARITY_RAD;
  }
  else if (fabs(fabs(e[1]) - M_PI_2) < tolerance)
  {
    e[0] = ANGLE_ON_SINGULARITY_RAD;
    e[2] = 2 * theta_minus + ANGLE_ON_SINGULARITY_RAD;
  }
  else // Safe
  {
    e[0] = theta_plus - theta_minus;
    e[2] = theta_plus + theta_minus;
  }

  if (not_proper)
  {
    e[1] -= M_PI_2;
    e[2] *= epsilon;
  }

  // Normalize to [-pi, pi]
  for (uint8_t i = 0; i < 3; i++)
  {
    if (e[i] < -M_PI)
    {
      e[i] += 2.0 * M_PI;
    }
    else if (e[i] > M_PI)
    {
      e[i] -= 2 * M_PI;
    }
  }
}

void euler_to_quat(const double e[3], const euler_seq_t es, double q[3])
{
  const double c[3] = {cos(0.5 * e[0]), cos(0.5 * e[1]), cos(0.5 * e[2])};
  const double s[3] = {sin(0.5 * e[0]), sin(0.5 * e[1]), sin(0.5 * e[2])};
  const double h[3] = {0.5 * e[0], 0.5 * e[1], 0.5 * e[2]};
  const double ch = cos(h[1]);
  const double sh = sin(h[1]);

  switch (es)
  {
  case EULER_XYZ:
  {
    q[0] = c[0] * c[1] * c[2] - s[0] * s[1] * s[2];
    q[1] = s[0] * c[1] * c[2] + c[0] * s[1] * s[2];
    q[2] = c[0] * s[1] * c[2] - s[0] * c[1] * s[2];
    q[3] = c[0] * c[1] * s[2] + s[0] * s[1] * c[2];
    break;
  }

  case EULER_XZY:
  {
    q[0] = c[0] * c[1] * c[2] + s[0] * s[1] * s[2];
    q[1] = s[0] * c[1] * c[2] - c[0] * s[1] * s[2];
    q[2] = c[0] * c[1] * s[2] - s[0] * s[1] * c[2];
    q[3] = c[0] * s[1] * c[2] + s[0] * c[1] * s[2];
    break;
  }

  case EULER_XYX:
  {
    q[0] = ch * cos(h[0] + h[2]);
    q[1] = ch * sin(h[0] + h[2]);
    q[2] = sh * cos(h[0] - h[2]);
    q[3] = sh * sin(h[0] - h[2]);
    break;
  }

  case EULER_XZX:
  {
    q[0] = ch * cos(h[0] + h[2]);
    q[1] = ch * sin(h[0] + h[2]);
    q[2] = sh * sin(-h[0] + h[2]);
    q[3] = sh * cos(-h[0] + h[2]);
    break;
  }

  case EULER_YXZ:
  {
    q[0] = c[0] * c[1] * c[2] + s[0] * s[1] * s[2];
    q[1] = c[0] * s[1] * c[2] + s[0] * c[1] * s[2];
    q[2] = s[0] * c[1] * c[2] - c[0] * s[1] * s[2];
    q[3] = c[0] * c[1] * s[2] - s[0] * s[1] * c[2];
    break;
  }

  case EULER_YZX:
  {
    q[0] = c[0] * c[1] * c[2] - s[0] * s[1] * s[2];
    q[1] = c[0] * c[1] * s[2] + s[0] * s[1] * c[2];
    q[2] = s[0] * c[1] * c[2] + c[0] * s[1] * s[2];
    q[3] = c[0] * s[1] * c[2] - s[0] * c[1] * s[2];
    break;
  }

  case EULER_YXY:
  {
    q[0] = ch * cos(h[0] + h[2]);
    q[1] = sh * cos(-h[0] + h[2]);
    q[2] = ch * sin(h[0] + h[2]);
    q[3] = sh * sin(-h[0] + h[2]);
    break;
  }

  case EULER_YZY:
  {
    q[0] = ch * cos(h[0] + h[2]);
    q[1] = sh * sin(h[0] - h[2]);
    q[2] = ch * sin(h[0] + h[2]);
    q[3] = sh * cos(h[0] - h[2]);
    break;
  }

  case EULER_ZXY:
  {
    q[0] = c[0] * c[1] * c[2] - s[0] * s[1] * s[2];
    q[1] = c[0] * s[1] * c[2] - s[0] * c[1] * s[2];
    q[2] = c[0] * c[1] * s[2] + s[0] * s[1] * c[2];
    q[3] = s[0] * c[1] * c[2] + c[0] * s[1] * s[2];
    break;
  }

  case EULER_ZYX:
  {
    q[0] = c[0] * c[1] * c[2] + s[0] * s[1] * s[2];
    q[1] = c[0] * c[1] * s[2] - s[0] * s[1] * c[2];
    q[2] = c[0] * s[1] * c[2] + s[0] * c[1] * s[2];
    q[3] = s[0] * c[1] * c[2] - c[0] * s[1] * s[2];

    break;
  }

  case EULER_ZXZ:
  {
    q[0] = ch * cos(h[0] + h[2]);
    q[1] = sh * cos(h[0] - h[2]);
    q[2] = sh * sin(h[0] - h[2]);
    q[3] = ch * sin(h[0] + h[2]);
    break;
  }

  case EULER_ZYZ:
  {
    q[0] = ch * cos(h[0] + h[2]);
    q[1] = sh * sin(-h[0] + h[2]);
    q[2] = sh * cos(-h[0] + h[2]);
    q[3] = ch * sin(h[0] + h[2]);
    break;
  }
  }
}

void euler_to_dcm(const double e[3], const euler_seq_t es, double m[3][3])
{
  const double c[3] = {cos(e[0]), cos(e[1]), cos(e[2])};
  const double s[3] = {sin(e[0]), sin(e[1]), sin(e[2])};

  switch (es)
  {
  case EULER_XYZ:
  {
    m[0][0] = c[1] * c[2];
    m[0][1] = c[2] * s[0] * s[1] + c[0] * s[2];
    m[0][2] = s[0] * s[2] - c[0] * c[2] * s[1];
    m[1][0] = -c[1] * s[2];
    m[1][1] = c[0] * c[2] - s[0] * s[1] * s[2];
    m[1][2] = c[2] * s[0] + c[0] * s[1] * s[2];
    m[2][0] = s[1];
    m[2][1] = -c[1] * s[0];
    m[2][2] = c[0] * c[1];
    break;
  }

  case EULER_XZY:
  {
    m[0][0] = c[1] * c[2];
    m[0][1] = c[0] * c[2] * s[1] + s[0] * s[2];
    m[0][2] = c[2] * s[0] * s[1] - c[0] * s[2];
    m[1][0] = -s[1];
    m[1][1] = c[0] * c[1];
    m[1][2] = c[1] * s[0];
    m[2][0] = c[1] * s[2];
    m[2][1] = -c[2] * s[0] + c[0] * s[1] * s[2];
    m[2][2] = c[0] * c[2] + s[0] * s[1] * s[2];
    break;
  }

  case EULER_XYX:
  {
    m[0][0] = c[1];
    m[0][1] = s[0] * s[1];
    m[0][2] = -c[0] * s[1];
    m[1][0] = s[1] * s[2];
    m[1][1] = c[0] * c[2] - c[1] * s[0] * s[2];
    m[1][2] = c[2] * s[0] + c[0] * c[1] * s[2];
    m[2][0] = c[2] * s[1];
    m[2][1] = -c[1] * c[2] * s[0] - c[0] * s[2];
    m[2][2] = c[0] * c[1] * c[2] - s[0] * s[2];
    break;
  }

  case EULER_XZX:
  {
    m[0][0] = c[1];
    m[0][1] = c[0] * s[1];
    m[0][2] = s[0] * s[1];
    m[1][0] = -c[2] * s[1];
    m[1][1] = c[0] * c[1] * c[2] - s[0] * s[2];
    m[1][2] = c[1] * c[2] * s[0] + c[0] * s[2];
    m[2][0] = s[1] * s[2];
    m[2][1] = -c[2] * s[0] - c[0] * c[1] * s[2];
    m[2][2] = c[0] * c[2] - c[1] * s[0] * s[2];
    break;
  }

  case EULER_YXZ:
  {
    m[0][0] = c[0] * c[2] + s[0] * s[1] * s[2];
    m[0][1] = c[1] * s[2];
    m[0][2] = -c[2] * s[0] + c[0] * s[1] * s[2];
    m[1][0] = c[2] * s[0] * s[1] - c[0] * s[2];
    m[1][1] = c[1] * c[2];
    m[1][2] = c[0] * c[2] * s[1] + s[0] * s[2];
    m[2][0] = c[1] * s[0];
    m[2][1] = -s[1];
    m[2][2] = c[0] * c[1];
    break;
  }

  case EULER_YZX:
  {
    m[0][0] = c[0] * c[1];
    m[0][1] = s[1];
    m[0][2] = -c[1] * s[0];
    m[1][0] = -c[0] * c[2] * s[1] + s[0] * s[2];
    m[1][1] = c[1] * c[2];
    m[1][2] = c[2] * s[0] * s[1] + c[0] * s[2];
    m[2][0] = c[2] * s[0] + c[0] * s[1] * s[2];
    m[2][1] = -c[1] * s[2];
    m[2][2] = c[0] * c[2] - s[0] * s[1] * s[2];
    break;
  }

  case EULER_YXY:
  {
    m[0][0] = c[0] * c[2] - c[1] * s[0] * s[2];
    m[0][1] = s[1] * s[2];
    m[0][2] = -c[2] * s[0] - c[0] * c[1] * s[2];
    m[1][0] = s[0] * s[1];
    m[1][1] = c[1];
    m[1][2] = c[0] * s[1];
    m[2][0] = c[1] * c[2] * s[0] + c[0] * s[2];
    m[2][1] = -c[2] * s[1];
    m[2][2] = c[0] * c[1] * c[2] - s[0] * s[2];
    break;
  }

  case EULER_YZY:
  {
    m[0][0] = c[0] * c[1] * c[2] - s[0] * s[2];
    m[0][1] = c[2] * s[1];
    m[0][2] = -c[1] * c[2] * s[0] - c[0] * s[2];
    m[1][0] = -c[0] * s[1];
    m[1][1] = c[1];
    m[1][2] = s[0] * s[1];
    m[2][0] = c[2] * s[0] + c[0] * c[1] * s[2];
    m[2][1] = s[1] * s[2];
    m[2][2] = c[0] * c[2] - c[1] * s[0] * s[2];
    break;
  }

  case EULER_ZXY:
  {
    m[0][0] = c[0] * c[2] - s[0] * s[1] * s[2];
    m[0][1] = c[2] * s[0] + c[0] * s[1] * s[2];
    m[0][2] = -c[1] * s[2];
    m[1][0] = -c[1] * s[0];
    m[1][1] = c[0] * c[1];
    m[1][2] = s[1];
    m[2][0] = c[2] * s[0] * s[1] + c[0] * s[2];
    m[2][1] = s[0] * s[2] - c[0] * c[2] * s[1];
    m[2][2] = c[1] * c[2];
    break;
  }

  case EULER_ZYX:
  {
    m[0][0] = c[1] * c[0];
    m[0][1] = c[1] * s[0];
    m[0][2] = -s[1];
    m[1][0] = s[2] * s[1] * c[0] - c[2] * s[0];
    m[1][1] = s[2] * s[1] * s[0] + c[2] * c[0];
    m[1][2] = s[2] * c[1];
    m[2][0] = c[2] * s[1] * c[0] + s[2] * s[0];
    m[2][1] = c[2] * s[1] * s[0] - s[2] * c[0];
    m[2][2] = c[2] * c[1];
    break;
  }

  case EULER_ZXZ:
  {
    m[0][0] = c[2] * c[0] - s[2] * c[1] * s[0];
    m[0][1] = c[2] * s[0] + s[2] * c[1] * c[0];
    m[0][2] = s[2] * s[1];
    m[1][0] = -s[2] * c[0] - c[2] * c[1] * s[0];
    m[1][1] = -s[2] * s[0] + c[2] * c[1] * c[0];
    m[1][2] = c[2] * s[1];
    m[2][0] = s[1] * s[0];
    m[2][1] = -s[1] * c[0];
    m[2][2] = c[1];
    break;
  }

  case EULER_ZYZ:
  {
    m[0][0] = c[0] * c[1] * c[2] - s[0] * s[2];
    m[0][1] = c[1] * c[2] * s[0] + c[0] * s[2];
    m[0][2] = -c[2] * s[1];
    m[1][0] = -c[2] * s[0] - c[0] * c[1] * s[2];
    m[1][1] = c[0] * c[2] - c[1] * s[0] * s[2];
    m[1][2] = s[1] * s[2];
    m[2][0] = c[0] * s[1];
    m[2][1] = s[0] * s[1];
    m[2][2] = c[1];
    break;
  }
  }
}

void dcm_to_euler(const double m[3][3], double e[3], const euler_seq_t es)
{
  switch (es)
  {
  case EULER_XYZ:
  {
    e[0] = atan2(-m[2][1], m[2][2]);
    e[1] = asin(m[2][0]);
    e[2] = atan2(-m[1][0], m[0][0]);
    break;
  }

  case EULER_XZY:
  {
    e[0] = atan2(m[1][2], m[1][1]);
    e[1] = asin(-m[1][0]);
    e[2] = atan2(m[2][0], m[0][0]);
    break;
  }

  case EULER_XYX:
  {
    e[0] = atan2(m[0][1], -m[0][2]);
    e[1] = acos(m[0][0]);
    e[2] = atan2(m[1][0], m[2][0]);
    break;
  }

  case EULER_XZX:
  {
    e[0] = atan2(m[0][2], m[0][1]);
    e[1] = acos(m[0][0]);
    e[2] = atan2(m[2][0], -m[1][0]);
    break;
  }

  case EULER_YXZ:
  {
    e[0] = atan2(m[2][0], m[2][2]);
    e[1] = asin(-m[2][1]);
    e[2] = atan2(m[0][1], m[1][1]);
    break;
  }

  case EULER_YZX:
  {
    e[0] = atan2(-m[0][2], m[0][0]);
    e[1] = asin(m[0][1]);
    e[2] = atan2(-m[2][1], m[1][1]);
    break;
  }

  case EULER_YXY:
  {
    e[0] = atan2(m[1][0], m[1][2]);
    e[1] = acos(m[1][1]);
    e[2] = atan2(m[0][1], -m[2][1]);
    break;
  }

  case EULER_YZY:
  {
    e[0] = atan2(m[1][2], -m[1][0]);
    e[1] = acos(m[1][1]);
    e[2] = atan2(m[2][1], m[0][1]);
    break;
  }

  case EULER_ZXY:
  {
    e[0] = atan2(-m[1][0], m[1][1]);
    e[1] = asin(m[1][2]);
    e[2] = atan2(-m[0][2], m[2][2]);
    break;
  }

  case EULER_ZYX:
  {
    e[0] = atan2(m[0][1], m[0][0]);
    e[1] = asin(-m[0][2]);
    e[2] = atan2(m[1][2], m[2][2]);
    break;
  }

  case EULER_ZXZ:
  {
    e[0] = atan2(m[2][0], -m[2][1]);
    e[1] = acos(m[2][2]);
    e[2] = atan2(m[0][2], m[1][2]);
    break;
  }

  case EULER_ZYZ:
  {
    e[0] = atan2(m[2][1], m[2][0]);
    e[1] = acos(m[2][2]);
    e[2] = atan2(m[1][2], -m[0][2]);
    break;
  }
  }
}
