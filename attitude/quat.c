#include <math.h>
#include <inttypes.h>

#include "quat.h"

#define PI_2 1.57079632679489661923
#define PI 2 * 1.57079632679489661923

#define QUAT_ANGLE_ON_SINGULARITY_RAD 0.0f
#define QUAT_COMPARISON_EPSILON 1e-5f
#define QUAT_DIVISION_EPSILON 1e-6f
#define EULER_SINGULARITY_EPSILON 1e-7f

// Initializes input array to the unit quaternion.
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

// Negates all the elements of the input quaternion.
void quat_neg(double q[4])
{
  q[0] = -q[0];
  q[1] = -q[1];
  q[2] = -q[2];
  q[3] = -q[3];
}

// Normalizes the input quaternion.
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

// Computes norm of the input quaternion.
double quat_norm(const double q[4])
{
  return sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
}

// Computes conjugate of the input quaternion.
void quat_conj(const double q[4], double q_conj[4])
{
  q_conj[0] = q[0];
  q_conj[1] = -q[1];
  q_conj[2] = -q[2];
  q_conj[3] = -q[3];
}

// Copies the contents of `q_src` to `q_dest`.
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

// Computes the dot product of the input quaternions.
double quat_dot(const double q1[4], const double q2[4])
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

// Quaternion product: q = q1 * q2
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

/**
 * @brief Computes the angle-axis representation corresponding to the input quaternion, where the
 *        rotation occurs by an angle `psi` [rad] about the unit axis vector `v`.
 *
 * @param q Input quaternion
 * @param psi Rotation angle in radians
 * @param v Unit vector representing the axis of rotation
 *
 * @return Status of the conversion; `true` if successful, `false` if there are infinite solutions,
 *         in which case the axis is chosen to be [1, 0, 0]'.
 */
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

/**
 * @brief Computes the rate of change of the quaternion using the input quaternion and angular velocity.
 *
 * The angular velocity vector `w` is assumed to be expressed in the body frame.
 *
 * @param q Input quaternion
 * @param w Angular velocity vector in the body frame
 * @param q_dot Output rate of change of the quaternion
 */
void quat_rate(const double q[4], const double w[3], double q_dot[4])
{
  q_dot[0] = 0.5 * (-w[0] * q[1] - w[1] * q[2] - w[2] * q[3]);
  q_dot[1] = 0.5 * (w[0] * q[0] + w[2] * q[2] - w[1] * q[3]);
  q_dot[2] = 0.5 * (w[1] * q[0] - w[2] * q[1] + w[0] * q[3]);
  q_dot[3] = 0.5 * (w[2] * q[0] + w[1] * q[1] - w[0] * q[2]);
}

// Rotates the input vector `v` using the quaternion `q` to produce the output vector `v_out`.
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

// Computes the DCM corresponding to the input quaternion.
void quat_to_dcm(const double q[4], double m[3][3])
{
  m[0][0] = q[0] * q[0] + q[1] * q[1] - q[2] * q[2] - q[3] * q[3];
  m[0][1] = 2.0 * (q[1] * q[2] + q[0] * q[3]);
  m[0][2] = 2.0 * (q[1] * q[3] - q[0] * q[2]);
  m[1][0] = 2.0 * (q[1] * q[2] - q[0] * q[3]);
  m[1][1] = q[0] * q[0] - q[1] * q[1] + q[2] * q[2] - q[3] * q[3];
  m[1][2] = 2.0 * (q[2] * q[3] + q[0] * q[1]);
  m[2][0] = 2.0 * (q[1] * q[3] + q[0] * q[2]);
  m[2][1] = 2.0 * (q[2] * q[3] - q[0] * q[1]);
  m[2][2] = q[0] * q[0] - q[1] * q[1] - q[2] * q[2] + q[3] * q[3];
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

static inline double det_3x3(double m[3][3])
{
  return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
         m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
         m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}

// Compute greatest eigenvalue of traceless symmetric 4x4 matrix
double get_lambda_max(const double K[4][4], const int n)
{
  // Cache unique elements of the symmetric matrix
  const double K00 = K[0][0], K01 = K[0][1], K02 = K[0][2], K03 = K[0][3];
  const double K11 = K[1][1], K12 = K[1][2], K13 = K[1][3];
  const double K22 = K[2][2], K23 = K[2][3];
  const double K33 = K[3][3];

  // Cache frequently used square terms
  const double K03_sq = K03 * K03;
  const double K13_sq = K13 * K13;
  const double K23_sq = K23 * K23;
  const double K12_sq = K12 * K12;
  const double K01_sq = K01 * K01;
  const double K02_sq = K02 * K02;

  // Trace and adjugate terms. Eqn.(7)
  const double trB = K33;
  const double tradB = // trace(adj(B + tr(B)))
      (K11 + trB) * (K22 + trB) - K12_sq +
      (K00 + trB) * (K22 + trB) - K02_sq +
      (K00 + trB) * (K11 + trB) - K01_sq;
  const double b = -2.0f * trB * trB + tradB - (K03_sq + K13_sq + K23_sq);

  // Determinant of symmetric 4x4 matrix
  const double d =
      K03_sq * K12_sq - K00 * K33 * K12_sq - 2.0 * K01 * K03 * K23 * K12 +
      2.0 * K01 * K02 * K33 * K12 + 2.0 * K00 * K23 * K12 * K13 -
      2.0 * K02 * K03 * K12 * K13 - K11 * K00 * K23_sq + K01_sq * K23_sq +
      2.0 * K11 * K02 * K03 * K23 - K11 * K02_sq * K33 + K02_sq * K13_sq -
      2.0 * K01 * K02 * K23 * K13 + K11 * K00 * K33 * K22 - K00 * K13_sq * K22 -
      K01_sq * K33 * K22 + 2.0 * K01 * K03 * K13 * K22 - K11 * K03_sq * K22;

  // Two quaternions
  if (n == 2)
  {
    const double sqrt_d = sqrt(d);
    const double g3 = sqrt(2 * sqrt_d - b);  // Eqn.(15)
    const double g4 = sqrt(-2 * sqrt_d - b); // Eqn.(15)
    return (g3 + g4) / 2.0f;                 // Eqn.(14)
  }

  // trace(adj(K))
  double m1[3][3] = {{K11, K12, K13}, {K12, K22, K23}, {K13, K23, K33}};
  double m2[3][3] = {{K00, K02, K03}, {K02, K22, K23}, {K03, K23, K33}};
  double m3[3][3] = {{K00, K01, K03}, {K01, K11, K13}, {K03, K13, K33}};
  double m4[3][3] = {{K00, K01, K02}, {K01, K11, K12}, {K02, K12, K22}};
  const double tradK = det_3x3(m1) + det_3x3(m2) + det_3x3(m3) + det_3x3(m4);
  const double c = -tradK; // Eqn.(7)

  // Eqn.(10)
  const double p = b * b / 9.0 + 4.0 * d / 3.0;
  const double q = b * b * b / 27.0 - 4.0 * d * b / 3.0 + c * c / 2.0;

  // Eigenvalue solution
  const double arg = q / pow(p, 1.5);
  const double theta = acos(fmax(-1.0, fmin(1.0, arg)));
  const double u1 = 2.0 * sqrt(p) * cos(theta / 3.0) + b / 3.0; // Eqn.(9)
  const double g1 = sqrt(u1 - b);                               // Eqn.(12)
  const double g2 = -2.0 * sqrt(u1 * u1 - 4.0 * d);             // Eqn.(12)
  return 0.5f * (g1 + sqrt(-u1 - b - g2));                      // Eqn.(13)
}

/**
 * @brief Computes the quaternion corresponding to the largest eigenvalue of M.
 *
 * @param[in] M 4x4 symmetric matrix to be eigenanalyzed.
 * @param[in] w_sum Sum of the positive weights for each quaternion.
 * @param[in] n Number of quaternions.
 * @param[out] q Output quaternion (unit eigenvector).
 *
 * @note M must first be converted to a traceless matrix K for ESOQ2 to work correctly.
 */
void esoq2(const double M[4][4], const double w_sum, const double n, double q[4])
{
  double K[4][4] = {0};

  // Compute K because ESOQ2 requires a traceless matrix. Eqn.(11)
  for (int i = 0; i < 4; i++)
  {
    for (int j = 0; j < 4; j++)
    {
      K[i][j] = 4.0 * M[i][j];
    }
  }

  K[0][0] -= w_sum;
  K[1][1] -= w_sum;
  K[2][2] -= w_sum;
  K[3][3] -= w_sum;

  const double lambda = get_lambda_max(K, n);

  // Extract components from K matrix
  const double S00 = K[0][0] - lambda;
  const double S11 = K[1][1] - lambda;
  const double S22 = K[2][2] - lambda;
  const double S01 = K[0][1];
  const double S02 = K[0][2];
  const double S12 = K[1][2];

  const double z0 = K[0][3];
  const double z1 = K[1][3];
  const double z2 = K[2][3];
  const double t = K[3][3] - lambda;

  // Entries of matrix M in Eqn.(20) computed using Eqn.(19)
  const double ma = S00 * t - z0 * z0;
  const double mb = S11 * t - z1 * z1;
  const double mc = S22 * t - z2 * z2;
  const double mx = S01 * t - z0 * z1;
  const double my = S02 * t - z0 * z2;
  const double mz = S12 * t - z1 * z2;

  double e[3] = {mb * mc - mz * mz, ma * mc - my * my, ma * mb - mx * mx};
  double e_abs[3] = {fabs(e[0]), fabs(e[1]), fabs(e[2])};

  // Choose optimal principal axis from three choices. Eqn.(21)
  if ((e_abs[0] > e_abs[1]) && (e_abs[0] > e_abs[2]))
  {
    e[1] = my * mz - mx * mc;
    e[2] = mx * mz - my * mb;
  }
  else if (e_abs[1] > e_abs[2])
  {
    e[0] = my * mz - mx * mc;
    e[2] = mx * my - mz * ma;
  }
  else
  {
    e[0] = mx * mz - my * mb;
    e[1] = mx * my - mz * ma;
  }

  // Normalize e
  const double norm_e_inv = 1.0 / sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2]);
  e[0] *= norm_e_inv;
  e[1] *= norm_e_inv;
  e[2] *= norm_e_inv;

  // Find maximum absolute value among the elements of z and t
  const double abs_z0 = fabs(z0);
  const double abs_z1 = fabs(z1);
  const double abs_z2 = fabs(z2);
  const double abs_t = fabs(t);

  double xk = 0.0;
  double yk = 0.0;

  if ((abs_z0 > abs_z1) && (abs_z0 > abs_z2) && (abs_z0 > abs_t))
  {
    xk = z0;
    yk = S00 * e[0] + S01 * e[1] + S02 * e[2];
  }
  else if ((abs_z1 > abs_z2) && (abs_z1 > abs_t))
  {
    xk = z1;
    yk = S01 * e[0] + S11 * e[1] + S12 * e[2];
  }
  else if (abs_z2 > abs_t)
  {
    xk = z2;
    yk = S02 * e[0] + S12 * e[1] + S22 * e[2];
  }
  else
  {
    xk = t;
    yk = z0 * e[0] + z1 * e[1] + z2 * e[2];
  }

  // Eqn.(28)
  const double h = sqrt(xk * xk + yk * yk);
  const double h_inv = 1.0 / h;
  const double sph = xk * h_inv;
  const double cph = -yk * h_inv;

  // I have no idea why, but rearranging seems to work.
  q[0] = sph * e[0];
  q[1] = sph * e[1];
  q[2] = sph * e[2];
  q[3] = cph;
}

// Weighted average n quaternions as described in Ref.[7].
void quat_mean(const double *q[], const int n, const double *w, double qm[4])
{
  if (n == 2)
  {
    const double *q1 = q[0];
    const double *q2 = q[1];
    const double w1 = w[0];
    const double w2 = w[1];

    double dot12 = quat_dot(q1, q2);
    double z = sqrt((w1 - w2) * (w1 - w2) + 4.0 * w1 * w2 * dot12 * dot12); // Eqn.(18)

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

    return;
  }

  double M[4][4] = {{0.0}};
  double w_sum = 0.0;

  for (int i = 0; i < n; i++)
  {
    w_sum += w[i];
  }

  // Compute the symmetric matrix M.  Eqn.(12)
  for (int i = 0; i < n; i++)
  {
    double qi[4] = {q[i][0], q[i][1], q[i][2], q[i][3]};

    for (int j = 0; j < 4; j++)
    {
      for (int k = 0; k < 4; k++)
      {
        M[j][k] += w[i] * qi[j] * qi[k];
      }
    }
  }

  esoq2(M, w_sum, n, qm);
  quat_canonize(qm);
}

/**
 * @brief Discrete-time propagation of quaternion kinematics.
 *
 * Solves the quaternion differential equation as expressed in `quat_rate()` using a discrete-time
 * approximation.
 *
 * @param q      Input quaternion (unit quaternion).
 * @param w      Angular velocity vector [rad/s].
 * @param dt     Time step [s].
 * @param q_out  Output quaternion after propagation.
 *
 * Reference:
 *   Crassidis, Junkins - Optimal estimation of dynamic systems (2011)
 */
void quat_prop(const double q[4], const double w[3], double dt, double q_out[4])
{
  double n = sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);

  if (n > QUAT_DIVISION_EPSILON)
  {
    double half_theta = 0.5 * n * dt;
    double c = cos(half_theta);
    double s = sin(half_theta) / n;

    double x = w[0] * s;
    double y = w[1] * s;
    double z = w[2] * s;

    q_out[0] =  c * q[0] + z * q[1] - y * q[2] + x * q[3];
    q_out[1] = -z * q[0] + c * q[1] + x * q[2] + y * q[3];
    q_out[2] =  y * q[0] - x * q[1] + c * q[2] + z * q[3];
    q_out[3] = -x * q[0] - y * q[1] - z * q[2] + c * q[3];
  }
  else
  {
    // Assume no rotation
    for (int i = 0; i < 4; ++i)
    {
      q_out[i] = q[i];
    }
  }

  quat_normalize(q_out);
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
  if (quat_equal(q1, q2))
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
  dot = (dot < -1.0) ? -1.0 : (dot > 1.0) ? 1.0
                                          : dot;

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
