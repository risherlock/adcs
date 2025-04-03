#include <math.h>
#include <inttypes.h>

#include "quat.h"

#define PI_2 1.57079632679489661923
#define PI 2 * 1.57079632679489661923

#define QUAT_ANGLE_ON_SINGULARITY_RAD 0.0f
#define QUAT_COMPARISON_EPSILON 1e-5f
#define QUAT_DIVISION_EPSILON 1e-9f

void quat_unit(double q[4])
{
  q[0] = 1.0;
  q[1] = 0.0;
  q[2] = 0.0;
  q[3] = 0.0;
}

void quat_normalize(double q[4])
{
  double norm = quat_norm(q);

  if (norm > QUAT_DIVISION_EPSILON)
  {
    double inv_norm = 1.0 / norm;
    q[0] = q[0] * inv_norm;
    q[1] = q[1] * inv_norm;
    q[2] = q[2] * inv_norm;
    q[3] = q[3] * inv_norm;
  }
}

double quat_norm(const double q[4])
{
  return sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
}

void quat_conj(const double q[4], double q_conj[4])
{
  q_conj[0] = q[0];
  q_conj[1] = -q[1];
  q_conj[2] = -q[2];
  q_conj[3] = -q[3];
}

void quat_copy(const double q_src[4], double q_dest[4])
{
  q_dest[0] = q_src[0];
  q_dest[1] = q_src[1];
  q_dest[2] = q_src[2];
  q_dest[3] = q_src[3];
}

bool quat_equal(const double q1[4], const double q2[4])
{
  bool q0_eq = fabs(q1[0] - q2[0]) <= QUAT_COMPARISON_EPSILON;
  bool q1_eq = fabs(q1[1] - q2[1]) <= QUAT_COMPARISON_EPSILON;
  bool q2_eq = fabs(q1[2] - q2[2]) <= QUAT_COMPARISON_EPSILON;
  bool q3_eq = fabs(q1[3] - q2[3]) <= QUAT_COMPARISON_EPSILON;

  return q0_eq && q1_eq && q2_eq && q3_eq;
}

double quat_dot(const double q1[4], const double q2[2])
{
  return q1[0] * q2[0] + q1[1] * q2[1] + q1[2] * q2[2] + q1[3] * q2[3];
}

// q = q1 * q2
void quat_prod(const double q1[4], const double q2[4], double q[4])
{
  q[0] = q1[0] * q2[0] - q1[1] * q2[1] - q1[2] * q2[2] - q1[3] * q2[3];
  q[1] = q1[0] * q2[1] + q1[1] * q2[0] + q1[2] * q2[3] - q1[3] * q2[2];
  q[2] = q1[0] * q2[2] - q1[1] * q2[3] + q1[2] * q2[0] + q1[3] * q2[1];
  q[3] = q1[0] * q2[3] + q1[1] * q2[2] - q1[2] * q2[1] + q1[3] * q2[1];
}

// qe = q1 - q2
void quat_err(const double q1[4], const double q2[4], double qe[4])
{
  double q_conj[4];
  quat_conj(q1, q_conj);
  quat_prod(q_conj, q2, qe);
}

bool quat_to_axan(const double q[4], double *psi, double v[3])
{
  bool status = true;
  *psi = 2.0 * acos(q[0]);
  double divisor = sqrt(1.0 - q[0] * q[0]);

  if (divisor > QUAT_DIVISION_EPSILON)
  {
    double inv_divisor = 1.0 / divisor;
    v[0] = q[1] * inv_divisor;
    v[1] = q[2] * inv_divisor;
    v[2] = q[3] * inv_divisor;
  }
  else // Infinite solutions
  {
    // Arbitrary axis
    v[0] = 1;
    v[1] = 0;
    v[2] = 0;

    status = false;
  }

  return status;
}

void quat_rotate(const double q[4], const double v[3], double v_out[3])
{
  double q0_sq = q[0] * q[0];
  double q1_sq = q[1] * q[1];
  double q2_sq = q[2] * q[2];
  double q3_sq = q[3] * q[3];
  double q0q1x2 = 2.0 * q[0] * q[1];
  double q0q2x2 = 2.0 * q[0] * q[2];
  double q0q3x2 = 2.0 * q[0] * q[3];
  double q1q2x2 = 2.0 * q[1] * q[2];
  double q1q3x2 = 2.0 * q[1] * q[3];
  double q2q3x2 = 2.0 * q[2] * q[3];

  v_out[0] = (q0_sq + q1_sq - q2_sq - q3_sq) * v[0] + (q1q2x2 - q0q3x2) * v[1] + (q1q3x2 + q0q2x2) * v[2];
  v_out[1] = (q1q2x2 + q0q3x2) * v[0] + (q0_sq - q1_sq + q2_sq - q3_sq) * v[1] + (q2q3x2 - q0q1x2) * v[2];
  v_out[2] = (q1q3x2 - q0q2x2) * v[0] + (q2q3x2 + q0q1x2) * v[1] + (q0_sq - q1_sq - q2_sq + q3_sq) * v[2];
}

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
    e[0] = QUAT_ANGLE_ON_SINGULARITY_RAD;
    e[2] = 2 * theta_plus - QUAT_ANGLE_ON_SINGULARITY_RAD;
  }
  else if (fabs(fabs(e[1]) - PI_2) < tolerance)
  {
    e[0] = QUAT_ANGLE_ON_SINGULARITY_RAD;
    e[2] = 2 * theta_minus + QUAT_ANGLE_ON_SINGULARITY_RAD;
  }
  else // Safe
  {
    e[0] = theta_plus - theta_minus;
    e[2] = theta_plus + theta_minus;
  }

  if (not_proper)
  {
    e[1] -= PI_2;
    e[2] *= epsilon;
  }

  // Normalize to [-pi, pi]
  for (uint8_t i = 0; i < 3; i++)
  {
    if (e[i] < -PI)
    {
      e[i] += 2.0 * PI;
    }
    else if (e[i] > PI)
    {
      e[i] -= 2 * PI;
    }
  }
}
