// Basic quaternion manipulations.
// 2024-12-22

#ifndef _ATTITUDE_QUAT_H_
#define _ATTITUDE_QUAT_H_

#include <stdbool.h>
#include "euler.h"

void quat_unit(double q[4]);
void quat_normalize(double q[4]);
double quat_norm(const double q[4]);
void quat_to_dcm(const double q[4], double m[3][3]);
void quat_conj(const double q[4], double q_conj[4]);
void quat_copy(const double q_src[4], double q_dest[4]);
bool quat_equal(const double q1[4], const double q2[4]);
double quat_dot(const double q1[4], const double q2[2]);
bool quat_to_axan(const double q[4], double *psi, double v[3]);
void quat_prod(const double q1[4], const double q2[4], double q[4]);
void quat_err(const double q1[4], const double q2[4], double qe[4]);
void quat_rotate(const double q[4], const double v[3], double v_out[3]);
void quat_to_euler(const double q[4], double e[3], const euler_seq_t es);

/*
  Todos
  -----
  1. quat_mean(const double qn[][4], const int n, double q[4]);
     const double *qs[] = {q0, q1, q2, q3}; -> quat_mean(qs, 4, q)
  2. quat_slerp(const double q1[4], const double q2[4], const double k, double q[4]);
  3. double quat_dist(const double q1[4], const double q2[4]);
  4. void quat_kine(const double q[4], const double w[3], double q_dot[4]);
*/

#endif // quat.h
