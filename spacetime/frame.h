// rms (2025-02-03)

#ifndef _ADCS_FRAME_H_
#define _ADCS_FRAME_H_

#include "frame.h"
#include "time.h"

void frame_eci_to_ecef_dcm(const utc_t t, double dcm[3][3]);

void frame_ecef_to_ned(const double v_ecef[3], const double llr[3], double v_ned[3]);
void frame_ned_to_ecef(const double v_ned[3], const double llr[3], double v_ecef[3]);

// void frame_eci_to_ned();
// void frame_ned_to_eci();

// void frame_ned_to_ecef_dcm(const double llr[3]);
// void frame_ecef_to_ned_dcm(const double llr[3]);

#endif // frame.h
