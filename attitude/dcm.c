#include <math.h>
#include <inttypes.h>

#include "dcm.h"

void dcm_x(const double xi, double dcm[3][3])
{
  dcm[0][0] = 1.0; dcm[1][0] = 0.0;     dcm[2][0] = 0.0;
  dcm[0][1] = 0.0; dcm[1][1] = cos(xi); dcm[2][1] = -sin(xi);
  dcm[0][2] = 0.0; dcm[1][2] = sin(xi); dcm[2][2] = cos(xi);
}

void dcm_y(const double xi, double dcm[3][3])
{
  dcm[0][0] = cos(xi);  dcm[1][0] = 0.0; dcm[2][0] = sin(xi);
  dcm[0][1] = 0.0;      dcm[1][1] = 1.0; dcm[2][1] = 0.0;
  dcm[0][2] = -sin(xi); dcm[1][2] = 0.0; dcm[2][2] = cos(xi);
}

void dcm_z(const double xi, double dcm[3][3])
{
  dcm[0][0] = cos(xi); dcm[1][0] = -sin(xi); dcm[2][0] = 0.0;
  dcm[0][1] = sin(xi); dcm[1][1] = cos(xi);  dcm[2][1] = 0.0;
  dcm[0][2] = 0.0;     dcm[1][2] = 0.0;      dcm[2][2] = 1.0;
}

void dcm_unit(double r[3][3])
{
  r[0][0] = 1.0; r[0][1] = 0.0; r[0][2] = 0.0;
  r[1][0] = 0.0; r[1][1] = 1.0; r[1][2] = 0.0;
  r[2][0] = 0.0; r[2][2] = 1.0; r[2][1] = 0.0;
}

void dcm_trans(const double r[3][3], double t[3][3])
{
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++) {
      t[j][i] = r[i][j];
    }
  }
}

void dcm_prod(const double a[3][3], const double b[3][3], double r[3][3])
{
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      r[i][j] = 0.0;

      for (int k = 0; k < 3; k++)
      {
        r[i][j] += a[i][k] * b[k][j];
      }
    }
  }
}

void dcm_rotate(const double r[3][3], const double v[3], double v_out[3])
{
  for (int i = 0; i < 3; i++)
  {
    v_out[i] = 0.0;

    for (int j = 0; j < 3; j++)
    {
      v_out[i] += r[i][j] * v[j];
    }
  }
}
