#include<malloc.h>
#include<stdio.h>
#include<stdlib.h>

#define _USE_MATH_DEFINES
#include<math.h>


/* Function to calculate intersection of two planes using dip and dip azimuth
Accepts input as follows: intersect(dip_az_1, dip_1, dip_az_2, dip_2, &output_azimuth, &output_plunge)*/
void intersect(double az1, double d1, double az2, double d2, double* out_az, double* out_dip){
    double p1az, p1pl, p2az, p2pl, p1cosa, p1cosb, p1cosc, p2cosa, p2cosb, p2cosc, theta, icosa, icosb, icosc, azimuth, plunge;
    
    // Calculate pole to  plane 1
    if(az1 >= 0.0 && az1 < 180.0){
        p1az = az1 + 180.0;
    } else if(az1 >= 180.0 && az1 < 360.0) {
        p1az = az1 - 180.0;
    } else {
        printf("Azimuth 1 not valid; exiting.");
        exit(0);
    }
    // Calculate pole to plane 2
    if(az2 >= 0.0 && az2 < 180.0){
        p2az = az2 + 180.0;
    } else if(az2 >= 180.0 && az2 < 360.0) {
        p2az = az2 - 180.0;
    } else {
        printf("Azimuth 2 not valid; exiting.");
        exit(0);
    }

    // Calculate plunge of pole to plane 1
    if(d1 <= 90.0 && d1 >= 0.0){
        p1pl = 90.0 - d1;
    } else {
        printf("Dip 1 not valid; exiting.");
        exit(0);
    }
    // Calculate plunge of pole to plane 2
    if(d2 <= 90.0 && d2 >= 0.0){
        p2pl = 90.0 - d2;
    } else {
        printf("Dip 2 not valid; exiting.");
        exit(0);
    }

    // Calculate cos alpha, beta, gamma for plane 1 (radians)
    p1cosa = sin(p1az * (M_PI / 180.0)) * sin((90 - p1pl) * (M_PI / 180.0));
    p1cosb = cos(p1az * (M_PI / 180.0)) * sin((90 - p1pl) * (M_PI / 180.0));
    p1cosc = cos((90 - p1pl) * (M_PI / 180.0));
    // Calculate cos alpha, beta, gamma for plane 2 (radians)
    p2cosa /*cn*/ = sin(p2az * (M_PI / 180.0)) * sin((90 - p2pl) * (M_PI / 180.0));
    p2cosb /*ce*/ = cos(p2az * (M_PI / 180.0)) * sin((90 - p2pl) * (M_PI / 180.0));
    p2cosc /*cd*/ = cos((90 - p2pl) * (M_PI / 180.0));

    // Calculate theta (in radians)
    theta = acos(p1cosa * p2cosa + p1cosb * p2cosb + p1cosc * p2cosc);

    // Calculate orientation in 3D
    if(theta != 0){
        icosa = (p1cosb * p2cosc - p1cosc * p2cosb)/sin(theta);
        icosb = -(p1cosa * p2cosc - p1cosc * p2cosa)/sin(theta);
        icosc = (p1cosa * p2cosb - p1cosb * p2cosa)/sin(theta);
    } else {
        icosa = 0;
        icosb = 0;
        icosc = 0;
    }

    // Convert to lower hemisphere
    if(icosc < 0){
        icosa = icosa * -1.0;
        icosb = icosb * -1.0;
        icosc = icosc * -1.0;
    }

    // Calculate intersection azimuth
    if(d1 == 90.0 && d2 == 90.0){
        azimuth = 0;
    } else if(icosa == 0.0 && icosb == 0.0){
        azimuth = 0.0;
    } else if(icosa < 0.0 && icosb >= 0.0){
        azimuth = 450.0 - (atan2(icosb,icosa) * 180.0/M_PI);
    } else {
        azimuth = 90.0 - (atan2(icosb,icosa) * 180.0/M_PI);
    }
    // Calculate intesection plunge
    if(d1 == 90.0 && d2 == 90.0){
        plunge = 90;
    } else {
        plunge = 90.0 - (acos(icosc) * 180.0/M_PI);
    }
    

    // Write results to pointers
    *out_az = azimuth;
    *out_dip = plunge;
}

double az, dip, az1, az2, d1, d2;
int i;

int main(int argc, char *argv[]){
    az1 = atof(argv[1]);
    az2 = atof(argv[3]);
    d1 = atof(argv[2]);
    d2 = atof(argv[4]);
    intersect(az1, d1, az2, d2, &az, &dip);
    printf("Azimuth: %.3f\tPlunge: %.3f", az, dip);
    return 0;
}