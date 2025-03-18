#include "igrf.h"
#include <math.h>

// World Geodetic System 1984 (WGS84) params
static const double wgs84_f = 1 / 298.257223563;
static const double wgs84_a = 6378.137;
static const double wgs84_b = wgs84_a * (1 - wgs84_f);

// Decimal years since January 1, IGRF_START_YEAR
double get_years(const igrf_time_t dt)
{
  if ((dt.year < IGRF_START_YEAR) || (dt.year >= IGRF_END_YEAR) ||
      (dt.month < 1) || (dt.month > 12) ||
      (dt.day < 1) || (dt.hour > 23) ||
      (dt.minute > 59) || (dt.second > 59))
  {
    return -1;
  }

  int days_in_month[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  int is_leap = (((dt.year % 4) == 0) && (((dt.year % 100) != 0) || ((dt.year % 400) == 0)));

  // Adjust for leap year
  if (is_leap)
  {
    days_in_month[1] = 29;
  }

  // Check valid day in the month
  if (dt.day > days_in_month[dt.month - 1])
  {
    return -1;
  }

  // Days since IGRF_START_YEAR
  int years = dt.year - IGRF_START_YEAR;
  int days_arr[] = {0, 31, 59, 90, 120, 151, 182, 212, 243, 273, 304, 334};
  double days = days_arr[dt.month - 1] + dt.day + (dt.month > 2 ? is_leap : 0);
  double hours = dt.hour + (dt.minute / 60.0) + (dt.second / 3600.0);
  int total_days = is_leap ? 366 : 365;

  // Decimal years
  return years + days / total_days + hours / 24.0;
}

/**
 * @brief Computes magnetic field strength [nT] in specified coordinate frame.
 *
 * @param t Time between January 1, IGRF_START_YEAR and December 31, IGRF_END_YEAR.
 * @param x Spatial coordinates: {latitude [deg], longitude [deg], altitude [km]}.
 *          - For geodetic frame, altitude represents height above Earth's surface.
 *          - For geocentric frame, altitude represents radius from the center of the Earth.
 * @param b Output magnetic field intensity in specified frame.
 * @param f The frame of reference for the input coordinates and the output field values.
 *
 * @return false if the time is out of range; true otherwise.
 */
bool igrf(const igrf_time_t t, const double x[3], const igrf_frame_t f, double b[3])
{
  const double a = 6371.2;                  // Radius of Earth [km]
  const double theta = M_PI_2 - x[0] * D2R; // Colattitude [rad]
  const double phi = x[1] * D2R;            // Longitude [rad]

  double cd, sd, r;
  double ct = cos(theta);
  double st = sin(theta);

  // Geodetic to geocentric conversion
  // https://github.com/wb-bgs/m_IGRF
  if (f == IGRF_GEODETIC)
  {
    // Radius
    const double h = x[2];
    const double rho = hypot(wgs84_a * st, wgs84_b * ct);
    r = sqrt(h * h + 2 * h * rho + (pow(wgs84_a, 4) * st * st + pow(wgs84_b, 4) * ct * ct) / (rho * rho));

    // Latitude
    cd = (h + rho) / r;
    sd = (wgs84_a * wgs84_a - wgs84_b * wgs84_b) / rho * ct * st / r;
    const double temp = ct;
    ct = cd * ct - sd * st;
    st = cd * st + sd * temp;
  }
  else if (f == IGRF_GEOCENTRIC)
  {
    r = x[2]; // Radius
  }

  // Avoid singularity on pole
  const double epsilon = 1e-8;

  if (st < epsilon && st > -epsilon)
  {
    st = epsilon;
  }

  double years = get_years(t);

  if (years < 0)
  {
    return false;
  }

  // [a] Re-occurring power factors
  // Optimizations [a] and [b] by Alar Leibak.
  double ar_pow[IGRF_DEGREE + 1];
  const double ar = a / r;
  ar_pow[0] = ar * ar * ar;

  for (uint8_t i = 1; i <= IGRF_DEGREE; i++)
  {
    ar_pow[i] = ar_pow[i - 1] * ar;
  }

  // [b] Re-occurring sines and cosines
  double sines[IGRF_DEGREE + 1], cosines[IGRF_DEGREE + 1];
  sines[0] = 0;
  cosines[0] = 1;
  sines[1] = sin(phi);
  cosines[1] = cos(phi);

  for (uint8_t i = 2; i <= IGRF_DEGREE; i++)
  {
    if (i & 1)
    {
      sines[i] = sines[i - 1] * cosines[1] + cosines[i - 1] * sines[1];
      cosines[i] = cosines[i - 1] * cosines[1] - sines[i - 1] * sines[1];
    }
    else // even
    {
      sines[i] = 2 * sines[i >> 1] * cosines[i >> 1];
      cosines[i] = 2 * cosines[i >> 1] * cosines[i >> 1] - 1;
    }
  }

  // Associated Legendre polynomials and its derivative
  double pnm = 0.0f, dpnm = 0.0f; // (n, m)
  double p11 = 1.0f, dp11 = 0.0f; // (n, n)
  double p10 = 1.0f, dp10 = 0.0f; // (n-1, m)
  double p20 = 0.0f, dp20 = 0.0f; // (n-2, m)

  // Field components: radial, theta, and phi
  double br = 0.0f, bt = 0.0f, bp = 0.0f;

  for (uint8_t m = 0; m <= IGRF_DEGREE; m++)
  {
    for (uint8_t n = 1; n <= IGRF_DEGREE; n++)
    {
      if (m <= n)
      {
        if (n == m)
        {
          pnm = st * p11;
          dpnm = st * dp11 + ct * p11;

          p11 = pnm; dp11 = dpnm;
          p20 = 0.0; dp20 = 0.0;
        }
        else
        {
          double Knm = 0.0;

          if (n > 1)
          {
            Knm = (pow(n - 1, 2) - pow(m, 2)) / ((2.0f * n - 1) * (2.0f * n - 3));
          }

          pnm = ct * p10 - Knm * p20;
          dpnm = ct * dp10 - st * p10 - Knm * dp20;
        }

        p20 = p10; dp20 = dp10;
        p10 = pnm; dp10 = dpnm;

        // Linear interpolation of g and h
        const int k = (0.5 * n * (n + 1) + m) - 1;
        const double g = g_val[k] + g_sv[k] * years;
        const double h = h_val[k] + h_sv[k] * years;

        if (m == 0)
        {
          const double temp = ar_pow[n - 1] * g;
          br += temp * (n + 1.0) * pnm;
          bt -= temp * dpnm;
        }
        else
        {
          const double hsin = h * sines[m];
          const double hcos = h * cosines[m];
          const double gsin = g * sines[m];
          const double gcos = g * cosines[m];
          const double temp = ar_pow[n - 1] * (gcos + hsin);

          br += temp * (n + 1.0) * pnm;
          bt -= temp * dpnm;
          bp -= ar_pow[n - 1] * m * (-gsin + hcos) * pnm;
        }
      }
    }
  }

  bp = bp / st;

  // Geocentric NED
  b[0] = -bt;
  b[1] = bp;
  b[2] = -br;

  // Geocentric to geodetic NED
  if (f == IGRF_GEODETIC)
  {
    double temp = b[0];
    b[0] = cd * b[0] + sd * b[2];
    b[2] = cd * b[2] - sd * temp;
  }

  return true;
}

// Magnetic inclination [rad]
double igrf_inc(const double b[3])
{
  const double h = sqrt(b[0] * b[0] + b[1] * b[1]);
  return atan(b[2] / h);
}

// Magnetic declination [rad]
double igrf_dec(const double b[3])
{
  return atan2(b[1], b[0]);
}

// Magnitude of magnetic field
double igrf_mag(const double b[3])
{
  return sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);
}
