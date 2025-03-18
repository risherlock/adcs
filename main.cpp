#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include "TLE.h"
#include "frame.h"
#include "attitude.h"
#include "time.h"

typedef struct VERIN {
    char line1[70];
    char line2[70];
    double startmin;
    double stepmin;
    double stopmin;
} VERIN;

void convert_position_to_latlong(double r[3], utc_t t, double lla[3]) {
    double dcm[3][3];
    // convert eci to ecef
    frame_eci_to_ecef_dcm(t,dcm);
    double v_ecef[3];
    dcm_rotate(dcm, r, v_ecef);

    const double ecef[3] = {v_ecef[0], v_ecef[1], v_ecef[2]};

    frame_ecef_to_lla(ecef, lla);
}

void compute_magnetic_field_at_position(double r[3], utc_t datetime, double b[3]) {
    double lla[3];
    convert_position_to_latlong(r, datetime,lla);
    // b=igrf(datetime, lla[0], lla[1], lla[2]);
}


void runVER(VERIN *verins)
{
    FILE *in_file = NULL;
    TLE tle;
    double r[3];
    double v[3];
    double rv[3];
    double vv[3];
    utc_t t;
    double b[3];
    double mins = 0;
    char line[256];
    std::ofstream out_file;
    int steps;

    steps = (verins->stopmin-verins->startmin)/verins->stepmin;

    out_file.open("out_file.csv");
    out_file << "timestamp," << "dis1,"<< "dis2," << "dis3,"<<"vel1,"<<"vel2,"<<"vel3"<<"b1"<<"b2"<<"b3";
    printf("steps: %d",steps);

    tle.parseLines(verins->line1,verins->line2);

    for(int i=0; i<=steps; i++) {
        mins = verins->startmin+i*verins->stepmin;
        tle.getRV(mins,r,v);
        compute_magnetic_field_at_position(r, tle.datetime, b);
    
        out_file<<"\n"<<mins<<","<<r[0]<<","<<r[1]<<","<<r[2]<<","<<v[0]<<","<<v[1]<<","<<v[2]<<","<<b[0]<<","<<b[1]<<","<<b[2];
    }


    if(out_file){
        out_file.close();
    }
}

int main(void) {
    VERIN verin; 

    // Manually populate the VERIN structure
    verin.stepmin = 10;
    verin.startmin = 0;
    verin.stopmin = 4320;
    strncpy(verin.line1, "1 25544U 98067A   20300.83097691  .00001534  00000-0  35580-4 0  9996", 69);
    verin.line1[69] = '\0'; 
    strncpy(verin.line2, "2 25544  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667", 69);
    verin.line2[69] = '\0';

    printf("File studied\n");
    
    runVER(&verin); // Pass the address of the stack-allocated structure
    printf("Verins ran\n");
    
    return 0;
}
