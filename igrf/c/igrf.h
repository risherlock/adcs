/**
 * @brief Implementation of International Geomagnetic Reference Field (IGRF-14) Model.
 * @cite Davis - Mathematical Modeling of Earth's Magnetic Field (2004)
 * @author risherlock
 * @date 2021-12-26
 */

#ifndef _IGRF_H_
#define _IGRF_H_

#ifdef __cplusplus
extern "C"
{
#endif

#include "igrf_2025.h"

#include <stdbool.h>
#include <inttypes.h>

#define R2D 57.2957795131
#define D2R 0.01745329251
#define M_PI_2 1.57079632679489661923

typedef struct
{
  uint16_t year;
  uint8_t month;
  uint8_t day;
  uint8_t hour;
  uint8_t minute;
  uint8_t second;
} igrf_time_t;

typedef enum
{
  IGRF_GEODETIC,
  IGRF_GEOCENTRIC
} igrf_frame_t;

double igrf_mag(const double b[3]);
double igrf_inc(const double b[3]);
double igrf_dec(const double b[3]);
bool igrf(const igrf_time_t t, const double x[3], const igrf_frame_t f, double b[3]);

#ifdef __cplusplus
}
#endif
#endif // igrf.h
