#include <math.h>
#include "time.h"

double fix(double x)
{
  if (x > 0)
  {
    return floor(x);
  }

  return ceil(x);
}

/**
 * @brief UTC to Julian date.
 *
 * Number of decimal days since noon on November 24, 4714 BCE in the proleptic
 * Gregorian calendar, or January 1, 4713 BCE in the proleptic Julian calendar.
 *
 * @param t UTC
 * @return Julian date
 */
double time_julian_date(const utc_t t)
{
  double j0 = 367 * t.year - fix(7.0 * (t.year + fix((t.month + 9.0) / 12.0)) / 4.0) + fix(275.0 * t.month / 9.0) + t.day + 1721013.5;
  double ut = t.hour + t.minute / 60.0 + t.second / 3600.0;

  return j0 + ut / 24.0;
}

double normalize_zero_to_360(double angle)
{
  while (angle < 0)
  {
    angle += 360.0;
  }

  while (angle >= 360.0)
  {
    angle -= 360.0;
  }

  return angle;
}

// UTC to Greenwich sidereal time.
double time_greenwich_sidereal(const utc_t t)
{
  const double j0 = 367 * t.year - fix(7.0 * (t.year + fix((t.month + 9.0) / 12.0)) / 4.0) + fix(275.0 * t.month / 9.0) + t.day + 1721013.5;
  const double j = (j0 - 2451545.0) / 36525.0;
  const double j_sq = j * j;
  const double g0 = normalize_zero_to_360(100.4606184 + 36000.77004 * j + 0.000387933 * j_sq - 2.583e-8 * j_sq * j);
  const double ut = t.hour + t.minute / 60.0 + t.second / 3600.0;

  return g0 + 360.98564724 * ut / 24.0;
}

/**
 * @brief UTC to local sidereal time [deg].
 * @param t Time in UTC
 * @param el East longitude [deg]
 * @return Local sidereal time [deg]
 */
double time_local_sidereal_deg(const utc_t t, const double elon)
{
  const double lst = time_greenwich_sidereal(t) + elon;
  return lst - 360.0 * fix(lst / 360.0);
}

// UTC to local sidereal time [hr].
double time_local_sidereal_hr(const utc_t t, const double elon)
{
  return time_local_sidereal_deg(t, elon) / 15.0;
}

int isLeapYear(int year)
{
  if (year % 4 != 0)
  {
    return 0;
  }
  else if (year % 100 != 0)
  {
    return 1;
  }
  else if (year % 400 != 0)
  {
    return 0;
  }
  else
  {
    return 1;
  }
}

double getDaysInMonth(int month, int year){
  double daysInMonth[] = {31,28,31,30,31,30,31,31,30,31,30,31};

  if (isLeapYear(year) == 1)
  {
    daysInMonth[1] = 29;
  }

  return daysInMonth[month-1];
}

void calculateMYD(const utc_t tin, double days, utc_t *tout){
  int month = tin.month;
  int year = tin.year; 
  double daysInMonth = getDaysInMonth(month, year);

  while (days > daysInMonth)
  {
    days -= daysInMonth;
    month++;
    if (month > 12)
    {
      month = 1;
      year++;
    }
    daysInMonth = getDaysInMonth(month, year);
  }

  tout->year = year;
  tout->month = month;
  tout->day = days;
}


double add_time_to_utc(const utc_t tin, const double seconds, utc_t *tout)
{
  double t = tin.hour * 3600.0 + tin.minute * 60.0 + tin.second + seconds;
  tout->second = fmod(t, 60.0);
  t = (t - tout->second) / 60.0;
  tout->minute = fmod(t, 60.0);
  t = (t - tout->minute) / 60.0;
  tout->hour = fmod(t, 24.0);
  t = (t - tout->hour) / 24.0;

  double days = t + tin.day;

  if(tin.day>27) {
    calculateMYD(tin, days, tout);
  } else {
    tout->day = t + tin.day;
    tout->month = tin.month;
    tout->year = tin.year;
  }

  return 0;
}
