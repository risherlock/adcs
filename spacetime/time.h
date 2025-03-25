/**
 * @brief Standard time conversions
 * @date 2024-12-24
 * @cite [1] Curtis - Orbital mechanics for engineering students (2020)
 */

#ifndef _ADCS_TIME_H_
#define _ADCS_TIME_H_

#include <inttypes.h>

typedef struct
{
  uint16_t year;
  uint8_t month;
  uint8_t day;
  uint8_t hour;
  uint8_t minute;
  uint8_t second;
}utc_t;

double time_julian_date(const utc_t t);
double time_greenwich_sidereal(const utc_t t);
double time_local_sidereal_hr(const utc_t t, const double elon);
double time_local_sidereal_deg(const utc_t t, const double elon);
double add_time_to_utc(const utc_t tin, const double seconds, utc_t *tout);

#endif // time.h
