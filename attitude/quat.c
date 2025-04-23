#include <math.h>
#include <inttypes.h>

#include "quat.h"

#define PI_2 1.57079632679489661923
#define PI 2 * 1.57079632679489661923

#define QUAT_ANGLE_ON_SINGULARITY_RAD 0.0f
#define QUAT_COMPARISON_EPSILON 1e-5f
#define QUAT_DIVISION_EPSILON 1e-6f
#define EULER_SINGULARITY_EPSILON 1e-7f

void quat_unit(double q[4])
{
  q[0] = 1.0;
  q[1] = 0.0;
  q[2] = 0.0;
  q[3] = 0.0;
}

// Enforce q[0] >= 0 i.e. the canonical form.
void quat_canonize(double q[4])
{
  if (q[0] < 0)
  {
    quat_neg(q);
  }
}

void quat_neg(double q[4])
{
  q[0] = -q[0];
  q[1] = -q[1];
  q[2] = -q[2];
  q[3] = -q[3];
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
  q_conj[0] =  q[0];
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

// Returns true if quaternions represent the same rotation (q1 == +/-q2).
bool quat_equal(const double q1[4], const double q2[4])
{
  const double dot = quat_dot(q1, q2);
  return fabs(1.0 - fabs(dot)) <= QUAT_COMPARISON_EPSILON;
}

double quat_dot(const double q1[4], const double q2[2])
{
  return q1[0] * q2[0] + q1[1] * q2[1] + q1[2] * q2[2] + q1[3] * q2[3];
}

// Angular distance [rad] between two quaternions.
double quat_dist(const double q1[4], const double q2[4])
{
  double qe[4];
  quat_err(q1, q2, qe); // Rotation from q2 to q1
  quat_normalize(qe);   // Confirm domain for acos
  quat_canonize(qe);    // Ensure the shortest path

  return 2.0 * acos(fabs(qe[0]));
}

// q = q1 * q2
void quat_prod(const double q1[4], const double q2[4], double q[4])
{
  q[0] = q1[0] * q2[0] - q1[1] * q2[1] - q1[2] * q2[2] - q1[3] * q2[3];
  q[1] = q1[0] * q2[1] + q1[1] * q2[0] + q1[2] * q2[3] - q1[3] * q2[2];
  q[2] = q1[0] * q2[2] - q1[1] * q2[3] + q1[2] * q2[0] + q1[3] * q2[1];
  q[3] = q1[0] * q2[3] + q1[1] * q2[2] - q1[2] * q2[1] + q1[3] * q2[0];
}

/**
 * @brief Computes the error quaternion representing the rotation from q2 to q1.
 *
 * @note In the setting of real numbers, x = x1 - x2 represents the offset that converts x2 to x1.
 *       Analogously the quaternion q (q1 - q2) represents the rotation which rotates q2 to q1,
 *       i.e. q1 = q * q2 => q = q1 * inv(q2).
 */
void quat_err(const double q1[4], const double q2[4], double q[4])
{
  double q2_inv[4];
  quat_conj(q2, q2_inv);
  quat_prod(q1, q2_inv, q);
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
  if (fabs(e[1]) < EULER_SINGULARITY_EPSILON)
  {
    e[0] = QUAT_ANGLE_ON_SINGULARITY_RAD;
    e[2] = 2 * theta_plus - QUAT_ANGLE_ON_SINGULARITY_RAD;
  }
  else if (fabs(fabs(e[1]) - PI_2) < EULER_SINGULARITY_EPSILON)
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

// Weighted average of two quaternions (n=2 case)
void quat_mean(const double *q[], const int n, const double *w, double qm[4])
{
  if (n != 2)
  {
    // TODO
    return;
  }

  const double *q1 = q[0];
  const double *q2 = q[1];
  const double w1 = w[0];
  const double w2 = w[1];

  double dot12 = quat_dot(q1, q2);
  double dot11 = quat_dot(q1, q1);
  double z = sqrt((w1 - w2) * (w1 - w2) + 4.0 * w1 * w2 * dot11 * dot11); // Eqn.(18)

  // Eqn.(19)
  double temp = z * (w1 + w2 + z);
  double a = sqrt(w1 * (w1 - w2 + z) / temp);
  double b = sqrt(w2 * (w2 - w1 + z) / temp);
  double sign = (dot12 >= 0) ? 1.0 : -1.0;

  for (int i = 0; i < 4; i++)
  {
    qm[i] = a * q1[i] + sign * b * q2[i];
  }

  quat_normalize(qm);
  quat_canonize(qm);
}

/**
 * @brief Spherical linear interpolation (SLERP) between two quaternions.
 *
 * @param q1 Start quaternion (must be a unit quaternion)
 * @param q2 End quaternion (must be a unit quaternion)
 * @param u Interpolation parameter in range [0,1]
 * @param q Interpolated output quaternion
 *
 * @note When q1 and q2 are antipodal (q1 = -q2), the interpolation plane is chosen arbitrarily due
 *       to infinite possible solutions.
 *
 *       The parameter `u` determines the interpolation between q1 and q2. u = 0 and u = 1 returns
 *       exactly q1 and q2 (or its antipodal equivalent -q2), respectively. Intermediate values
 *       (e.g. u = 0.6) results in a quaternion representing 60% of the rotation from q1 toward q2
 *       along the shortest path.
 */
void quat_slerp(const double q1[4], const double q2[4], const double u, double q[4])
{
  // If q1 and q2 represent same orientation.
  if(quat_equal(q1, q2))
  {
    quat_copy(q1, q);
    return;
  }

  double qb[4] = {q2[0], q2[1], q2[2], q2[3]};
  double dot = quat_dot(q1, q2);

  // Ensure shortest angular path
  if (dot < 0.0)
  {
    quat_neg(qb);
    dot = -dot;
  }

  // Ensure dot is in valid domain for acos i.e. [-1, 1].
  dot = (dot < -1.0) ? -1.0 : (dot > 1.0) ? 1.0 : dot;

  const double theta = acos(fabs(dot));
  const double st = sqrt(1.0 - dot * dot);
  double a, b;

  // Fall back to LERP if sin(theta) -> 0.
  if (st <= QUAT_DIVISION_EPSILON)
  {
    a = 1.0 - u;
    b = u;
  }
  else // SLERP
  {
    const double stinv = 1 / st;
    const double ut = u * theta;

    a = sin(theta - ut) * stinv;
    b = sin(ut) * stinv;
  }

  q[0] = a * q1[0] + b * qb[0];
  q[1] = a * q1[1] + b * qb[1];
  q[2] = a * q1[2] + b * qb[2];
  q[3] = a * q1[3] + b * qb[3];

  quat_normalize(q);
}
