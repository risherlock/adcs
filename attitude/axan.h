// Basic axis-angle transformations.
// rms (2025-04-02 1:19:22 AM)

#ifndef _ATTITUDE_AXAN_H_
#define _ATTITUDE_AXAN_H_

void axan_to_quat(const double psi, const double v[3], double q[4]);

#endif // axan.h
