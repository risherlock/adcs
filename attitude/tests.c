#include <stdio.h>

#include "attitude.h"

void quaternion_tests(void);

int main(void)
{
  printf("Hello world!\n");

  double e[3] = {0.1, 0.2, 0.2};
  double m[3][3];
  euler_to_dcm(e, EULER_XZY, m);

  double q[4] = {1.0, 2.0, 3.0, 4.0};
  quat_normalize(q);
  quat_to_euler(q, e, EULER_XZY);
  printf("%f, %f, %f, %f\n", q[0], q[1], q[2], q[3]);
  printf("%f, %f, %f\n", e[0], e[1], e[2]);

  return 0;
}

void quaternion_tests(void)
{
  double q1[4];
  quat_unit(q1);

  double q2[4];
  const double psi = 1.27;
  const double v[3] = {1.0, 0.0, 0.0};
  quat_axis_angle(psi, v, q2);

  double v_rot[3];
  quat_rotate(q2, v, v_rot);
}
