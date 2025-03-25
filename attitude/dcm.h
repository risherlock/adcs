// Basic Direction Cosine Matrix manipulations.
// 2024-12-23

#ifndef _ATTITUDE_DCM_H_
#define _ATTITUDE_DCM_H_

void dcm_unit(double m[3][3]);
void dcm_x(const double xi, double m[3][3]);
void dcm_y(const double xi, double m[3][3]);
void dcm_z(const double xi, double m[3][3]);
void dcm_trans(const double m[3][3], double t[3][3]);
void dcm_prod(const double a[3][3], const double b[3][3], double m[3][3]);
void dcm_rotate(const double m[3][3], const double v[3], double v_out[3]);

#endif // dcm.h
