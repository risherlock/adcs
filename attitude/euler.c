#include <math.h>
#include <stdbool.h>
#include <inttypes.h>

#include "euler.h"

void euler_to_quat(const double e[3], const euler_seq_t es, double q[4])
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
