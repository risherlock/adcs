#include <math.h>

#include "iau06.h"
#include "iau06_data.h"

#define PI 2 * 1.57079632679489661923
#define D2R 0.01745329251

void iau06_compute_nutation_v(const double t[5], double nv[14])
{
  double mMoon = 485868.249036 + 1717915923.2178 * t[0] + 31.8792 * t[1] + 0.051635 * t[2] - 0.00024470 * t[3];
  double mSun = 1287104.793048 + 129596581.0481 * t[0] - 0.5532 * t[1] + 0.000136 * t[2] - 0.00001149 * t[3];
  double umMoon = 335779.526232 + 1739527262.8478 * t[0] - 12.7512 * t[1] - 0.001037 * t[2] + 0.00000417 * t[3];
  double dSun = 1072260.703692 + 1602961601.2090 * t[0] - 6.3706 * t[1] + 0.006593 * t[2] - 0.00003169 * t[3];
  double omegaMoon = 450160.398036 - 6962890.5431 * t[0] + 7.4722 * t[1] + 0.007702 * t[2] - 0.00005939 * t[3];

  double lMercury = 4.402608842 + 2608.7903141574 * t[0];
  double lVenus = 3.176146697 + 1021.3285546211 * t[0];
  double lEarth = 1.753470314 + 628.3075849991 * t[0];
  double lMars = 6.203480913 + 334.06124267 * t[0];
  double lJupiter = 0.599546497 + 52.9690962641 * t[0];
  double lSaturn = 0.874016757 + 21.329910496 * t[0];
  double lUranus = 5.481293872 + 7.4781598567 * t[0];
  double lNeptune = 5.311886287 + 3.8133035638 * t[0];
  double pa = 0.02438175 * t[0] + 0.00000538691 * t[1];

  double arcsec_to_rad = D2R / 3600.0;
  nv[0] = fmod(mMoon * arcsec_to_rad, 2 * PI);
  nv[1] = fmod(mSun * arcsec_to_rad, 2 * PI);
  nv[2] = fmod(umMoon * arcsec_to_rad, 2 * PI);
  nv[3] = fmod(dSun * arcsec_to_rad, 2 * PI);
  nv[4] = fmod(omegaMoon * arcsec_to_rad, 2 * PI);
  nv[5] = fmod(lMercury, 2 * PI);
  nv[6] = fmod(lVenus, 2 * PI);
  nv[7] = fmod(lEarth, 2 * PI);
  nv[8] = fmod(lMars, 2 * PI);
  nv[9] = fmod(lJupiter, 2 * PI);
  nv[10] = fmod(lSaturn, 2 * PI);
  nv[11] = fmod(lUranus, 2 * PI);
  nv[12] = fmod(lNeptune, 2 * PI);
  nv[13] = fmod(pa, 2 * PI);
}

double iau06_compute_x(const double nv[14], const double t[5])
{
  double q;
  double x = 0.0;

  for (int i = 0; i < 1306; i++)
  {
    q = 0.0;
    for (int k = 0; k < 14; k++)
    {
      q += iau06_cip_x0[i][k + 3] * nv[k];
    }
    x += (iau06_cip_x0[i][1] * sin(q) + iau06_cip_x0[i][2] * cos(q));
  }

  for (int i = 0; i < 253; i++)
  {
    q = 0.0;
    for (int k = 0; k < 14; k++)
    {
      q += iau06_cip_x1[i][k + 3] * nv[k];
    }
    x += (iau06_cip_x1[i][1] * sin(q) + iau06_cip_x1[i][2] * cos(q)) * t[0];
  }

  for (int i = 0; i < 36; i++)
  {
    q = 0.0;
    for (int k = 0; k < 14; k++)
    {
      q += iau06_cip_x2[i][k + 3] * nv[k];
    }
    x += (iau06_cip_x2[i][1] * sin(q) + iau06_cip_x2[i][2] * cos(q)) * t[1];
  }

  for (int i = 0; i < 4; i++)
  {
    q = 0.0;
    for (int k = 0; k < 14; k++)
    {
      q += iau06_cip_x3[i][k + 3] * nv[k];
    }
    x += (iau06_cip_x3[i][1] * sin(q) + iau06_cip_x3[i][2] * cos(q)) * t[2];
  }

  q = 0.0;
  for (int k = 0; k < 14; k++)
  {
    q += iau06_cip_x4[k + 3] * nv[k];
  }
  x += (iau06_cip_x4[1] * sin(q) + iau06_cip_x4[2] * cos(q)) * t[3];

  return x;
}

double iau06_compute_y(const double nv[14], const double t[5])
{
  double q;
  double y = 0.0;

  for (int i = 0; i < 962; i++)
  {
    q = 0.0;
    for (int k = 0; k < 14; k++)
    {
      q += iau06_cip_y0[i][k + 3] * nv[k];
    }
    y += (iau06_cip_y0[i][1] * sin(q) + iau06_cip_y0[i][2] * cos(q));
  }

  for (int i = 0; i < 277; i++)
  {
    q = 0.0;
    for (int k = 0; k < 14; k++)
    {
      q += iau06_cip_y1[i][k + 3] * nv[k];
    }
    y += (iau06_cip_y1[i][1] * sin(q) + iau06_cip_y1[i][2] * cos(q)) * t[0];
  }

  for (int i = 0; i < 30; i++)
  {
    q = 0.0;
    for (int k = 0; k < 14; k++)
    {
      q += iau06_cip_y2[i][k + 3] * nv[k];
    }
    y += (iau06_cip_y2[i][1] * sin(q) + iau06_cip_y2[i][2] * cos(q)) * t[1];
  }

  for (int i = 0; i < 5; i++)
  {
    q = 0.0;
    for (int k = 0; k < 14; k++)
    {
      q += iau06_cip_y3[i][k + 3] * nv[k];
    }
    y += (iau06_cip_y3[i][1] * sin(q) + iau06_cip_y3[i][2] * cos(q)) * t[2];
  }

  q = 0.0;
  for (int k = 0; k < 14; k++)
  {
    q += iau06_cip_y4[k + 3] * nv[k];
  }
  y += (iau06_cip_y4[1] * sin(q) + iau06_cip_y4[2] * cos(q)) * t[3];

  return y;
}

double iau06_compute_s(const double nv[14], const double t[5])
{
  double q;
  double s = 0.0;

  for (int i = 0; i < 33; i++)
  {
    q = 0.0;
    for (int k = 0; k < 8; k++)
    {
      int nut_index = (k < 5) ? k : (k + 2);
      q += iau06_cip_s0[i][k + 3] * nv[nut_index];
    }
    s += (iau06_cip_s0[i][1] * sin(q) + iau06_cip_s0[i][2] * cos(q));
  }

  for (int i = 0; i < 3; i++)
  {
    q = 0.0;
    for (int k = 0; k < 8; k++)
    {
      int nut_index = (k < 5) ? k : (k + 2);
      q += iau06_cip_s1[i][k + 3] * nv[nut_index];
    }
    s += (iau06_cip_s1[i][1] * sin(q) + iau06_cip_s1[i][2] * cos(q)) * t[0];
  }

  for (int i = 0; i < 25; i++)
  {
    q = 0.0;
    for (int k = 0; k < 8; k++)
    {
      int nut_index = (k < 5) ? k : (k + 2);
      q += iau06_cip_s2[i][k + 3] * nv[nut_index];
    }
    s += (iau06_cip_s2[i][1] * sin(q) + iau06_cip_s2[i][2] * cos(q)) * t[1];
  }

  for (int i = 0; i < 4; i++)
  {
    q = 0.0;
    for (int k = 0; k < 8; k++)
    {
      int nut_index = (k < 5) ? k : (k + 2);
      q += iau06_cip_s3[i][k + 3] * nv[nut_index];
    }
    s += (iau06_cip_s3[i][1] * sin(q) + iau06_cip_s3[i][2] * cos(q)) * t[2];
  }

  q = 0.0;
  for (int k = 0; k < 8; k++)
  {
    int nut_index = (k < 5) ? k : (k + 2);
    q += iau06_cip_s4[k + 3] * nv[nut_index];
  }
  s += (iau06_cip_s4[1] * sin(q) + iau06_cip_s4[2] * cos(q)) * t[3];

  return s;
}

#include <stdio.h>

void iau06_get_xys(const double tt, double *x, double *y, double *s)
{
  double t[5] = {tt, pow(tt, 2), pow(tt, 3), pow(tt, 4), pow(tt, 5)};

  double nv[14];
  iau06_compute_nutation_v(t, nv);

  double x0 = -16617 + 2004191898 * t[0] - 429782.9 * t[1] - 198618.34 * t[2] + 7.578 * t[3] + 5.9285 * t[4];
  double y0 = -6951 - 25896 * t[0] - 22407274.7 * t[1] + 1900.59 * t[2] + 1112.526 * t[3] + 0.1358 * t[4];
  double s0 = 94 + 3808.65 * t[0] - 122.68 * t[1] - 72574.11 * t[2] + 27.98 * t[3] + 15.62 * t[4];

  double x_temp = x0 + iau06_compute_x(nv, t);
  double y_temp = y0 + iau06_compute_y(nv, t);
  double s_temp = s0 + iau06_compute_s(nv, t);

  printf("x: %.20f\n", x_temp);
  printf("y: %.20f\n", y_temp);
  printf("s: %.20f\n", s_temp);

  *x = (x_temp * 1e-6 / 3600.0) * D2R;
  *y = (y_temp * 1e-6 / 3600.0) * D2R;
  *s = (s_temp * 1e-6 / 3600.0) * D2R - x_temp * y_temp / 2.0;
}
