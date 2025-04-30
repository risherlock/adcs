#include <math.h>
#include "axan.h"

/**
 * @brief Computes the quaternion corresponding to a rotation by an angle `psi` [rad] about the unit
 *        axis vector `v`.
 *
 * @param psi Rotation angle in radians
 * @param v Unit vector representing the axis of rotation
 * @param q Output quaternion
 */
void axan_to_quat(const double psi, const double v[3], double q[4])
{
  double hp = psi / 2.0;
  double shp = sin(hp);

  q[0] = cos(hp);
  q[1] = shp * v[0];
  q[2] = shp * v[1];
  q[3] = shp * v[2];
}
