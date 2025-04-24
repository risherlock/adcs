#include <math.h>
#include <inttypes.h>

#include "dcm.h"

void dcm_x(const double xi, double m[3][3])
{
  m[0][0] = 1.0; m[1][0] = 0.0;     m[2][0] = 0.0;
  m[0][1] = 0.0; m[1][1] = cos(xi); m[2][1] = -sin(xi);
  m[0][2] = 0.0; m[1][2] = sin(xi); m[2][2] = cos(xi);
}

void dcm_y(const double xi, double m[3][3])
{
  m[0][0] = cos(xi);  m[1][0] = 0.0; m[2][0] = sin(xi);
  m[0][1] = 0.0;      m[1][1] = 1.0; m[2][1] = 0.0;
  m[0][2] = -sin(xi); m[1][2] = 0.0; m[2][2] = cos(xi);
}

void dcm_z(const double xi, double m[3][3])
{
  m[0][0] = cos(xi); m[1][0] = -sin(xi); m[2][0] = 0.0;
  m[0][1] = sin(xi); m[1][1] = cos(xi);  m[2][1] = 0.0;
  m[0][2] = 0.0;     m[1][2] = 0.0;      m[2][2] = 1.0;
}

void dcm_unit(double m[3][3])
{
  m[0][0] = 1.0; m[0][1] = 0.0; m[0][2] = 0.0;
  m[1][0] = 0.0; m[1][1] = 1.0; m[1][2] = 0.0;
  m[2][0] = 0.0; m[2][1] = 0.0; m[2][2] = 1.0;
}

void dcm_trans(const double m[3][3], double t[3][3])
{
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++) {
      t[j][i] = m[i][j];
    }
  }
}

void dcm_prod(const double a[3][3], const double b[3][3], double m[3][3])
{
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      m[i][j] = 0.0;

      for (int k = 0; k < 3; k++)
      {
        m[i][j] += a[i][k] * b[k][j];
      }
    }
  }
}

void dcm_rotate(const double m[3][3], const double v[3], double v_out[3])
{
  for (int i = 0; i < 3; i++)
  {
    v_out[i] = 0.0;

    for (int j = 0; j < 3; j++)
    {
      v_out[i] += m[i][j] * v[j];
    }
  }
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


/**
 * @brief Time derivative of DCM for given angular velocity expressed in the body frame.
 *
 * @note  This function assumes that the input DCM represents a rotation from the inertial frame to
 *        the body frame.
 *
 * @warning  If the DCM represents a rotation from the body frame to the inertial frame, the
 *           following equation applies: m_dot = m x [w]. In this case, the angular velocity is
 *           expressed in the inertial frame.
 *
 * @todo Explore how same function would be used for the case in the warning.
 */
void dcm_rate(const double w[3], const double m[3][3], double m_dot[3][3])
{
  // m_dot = -[w] x m where, [w] is skew-symmetric matrix of w.
  const double wx[3][3] = {{0.0, w[2],  -w[1]}, {-w[2], 0.0, w[0]}, {w[1], -w[0], 0.0}};
  dcm_prod(wx, m, m_dot);
}
