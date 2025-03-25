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
