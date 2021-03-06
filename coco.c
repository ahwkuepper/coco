/***************************************************************************
 *   Copyright (C) 2015 by Andreas H.W. Kuepper                            *
 *   ahwkuepper@gmail.com                                                  *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
/***************************************************************************
 *   Converts positions and velocities between Cartesian, equatorial and   *
 *   galactic coordinate systems.                                          *
 ***************************************************************************/
/***************************************************************************
 *   Compile using the command: cc -o coco coco.c -lm                      *
 ***************************************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include"coco.h"


int main (int argc, const char * argv[])
{

    double x[3] = {0.0, 0.0, 0.0};	// galactocentric cartesian coordinates [pc]
    double v[3] = {0.0, 0.0, 0.0};	// galactocentric cartesian rest-frame velocity [km/s]
	double dsun = 6000.0;//8330.0;//102800.0;//23500.0;//23426.700;//30000;	//distance from the sun [pc]
	double vrsun = 232.0;//72.3;//-58.7;  //heliocentric radial velocity [km/s]
	double vr = 0.0; //radial velocity in LSR frame corrected for solar motion [km/s]
	double l = 309.0;//181.0;//28.74;//202.31;//Pal4 0.852;//152.45;	//galactic longitude [deg]
	double b = 15.0;//-45.0;//42.19;//71.80;//Pal4 45.860;//37.44; //galactic lattitude [deg]
	double lcosb = 0.0; //galactic longitude corrected for lattitude factor [deg]
	double RA = 0.0;	//right ascension [deg]
	double DEC = 0.0;	//declination [deg]
	double mu_alpha = 0.0; //proper motion in RA direction [mas/yr]
	double mu_alphacosdelta = -5.1;//-2.259280; //-2.414470; //proper motion in RA direction corrected for declination factor [mas/yr]
	double mu_delta = -3.6;//-2.229050;//-2.391090; //proper motion in DEC direction [mas/yr]
	double mu = 0.51; //transverse motion [mas/yr]
	double PA = 205.1;//223.667780; //Position angle of transverse motion [deg]
	double mu_l = 0.0; //proper motion in l direction [mas/yr]
	double mu_lcosb = 0.0; //proper motion in l direction corrected for lattitude factor [mas/yr]
	double mu_b = 0.0; //proper motion in b direction [mas/yr]

	int coordtype = 2; //(1 = equatorial to galactic and cartesian, 2 = galactic to equatorial and cartesian, 3 = cartesian to equatorial and galactic)
	int vcoordtype = 2; //(1 = mu & position angle; 2 = mu_alpha & mu_delta or mu_alphacos(delta) & mu_delta; 3 = cartesian velocities)
	
	printf("\nINPUT:\n");	
	printf("X = %f\t Y = %f\t Z = %f [pc]\t Rgal = %f\n", x[0],x[1],x[2], sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]));
    printf("Vx = %f\t Vy = %f\t Vz = %f\t V = %f [km/s]\n", v[0],v[1],v[2], sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]));
    printf("dsun = %f [pc]\t vrsun = %f\t vr = %f [km/s]\n", dsun,vrsun,vr);
    printf("l = %f\t b = %f\t lcos(b) = %f [deg]\n", l,b,lcosb);
    printf("RA = %f\t DEC = %f [deg]\n", RA,DEC);
    printf("mu_alpha = %f\t mu_alphacos(delta) = %f\t mu_delta = %f\t mu = %f [mas/yr]\t PA = %f [deg]\tmu_l = %f\t mu_lcos(b) = %f\t mu_b = %f\n", mu_alpha,mu_alphacosdelta,mu_delta,mu,PA, mu_l,mu_lcosb,mu_b);

	convert(x, v, &dsun, &vrsun, &vr, &l, &b, &lcosb, &RA, &DEC, &mu_alpha, &mu_alphacosdelta, &mu_delta, &mu, &PA, &mu_l, &mu_b, &mu_lcosb, coordtype, vcoordtype, vLSR);
	
	printf("\nOUTPUT:\n");	
	printf("X = %f\t Y = %f\t Z = %f [pc]\t Rgal = %f\n", x[0],x[1],x[2], sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]));
    printf("Vx = %f\t Vy = %f\t Vz = %f\t V = %f [km/s]\n", v[0],v[1],v[2], sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]));
    printf("dsun = %f [pc]\t vrsun = %f\t vr = %f [km/s]\n", dsun,vrsun,vr);
    printf("l = %f\t b = %f\t lcos(b) = %f [deg]\n", l,b,lcosb);
    printf("RA = %f\t DEC = %f [deg]\n", RA,DEC);
    printf("mu_alpha = %f\t mu_alphacos(delta) = %f\t mu_delta = %f\t mu = %f [mas/yr]\t PA = %f [deg]\tmu_l = %f\t mu_lcos(b) = %f\t mu_b = %f\n", mu_alpha,mu_alphacosdelta,mu_delta,mu,PA, mu_l,mu_lcosb,mu_b);
	
    return 0;

}


void convert(double *xtemp, double *vtemp, double *dsuntemp, double *vrsuntemp, double *vrtemp, double *ltemp, double *btemp, double *lcosbtemp, double *RAtemp, double *DECtemp, double *mu_alphatemp, double *mu_alphacosdeltatemp, double *mu_deltatemp, double *mutemp, double *PAtemp, double *mu_ltemp, double *mu_btemp, double *mu_lcosbtemp, int coordtype, int vcoordtype, double vLSRtemp){
	
	double x,y,z;        //galactic coordinates [kpc]
	double xsun = -rgalsun/1000.0;//galactocentric distance of sun [kpc]
	double dsun;			//heliocentric radial distance [kpc]
	double dxy;             //heliocentric distance in xy-plane [kpc]
	double vx,vy,vz;        //cluster velocity [km/s]
	double vrx,vry,vrz;     //cluster radial velocity in 3d coordinates [km/s]
	double vtx,vty,vtz;     //cluster tangential velocity in 3d coordinates [km/s]
	double t,d,a;
	double T[3][3], A[3][3], B[3][3];
	double TI[3][3];
	double detT;
	double RArad, DECrad;
	double brad, lrad, lcosbrad;
    double mu_alphacosdelta, mu_delta, mu, PArad, mu_l, mu_b, mu_lcosb;
    double vr, vrsun;
	double RAENP = 0.0, DECENP = PI/2.0;  //equatorial coordinates of equatorial north pole
	double xENP, yENP, zENP, dxyENP; //cartesian vector pointing to the equatorial north pole  
	double bENP, lENP; //galactic coordinates of ENP
	double FAK;
	double xdelta, ydelta, zdelta;
	double nx, ny, nz, dvt;
	double vrLSR, vrGSR;
	double C1, C2;
	
	//transformation matrix equ. -> gal. from Johnson & Soderblom (1987)
	t = PAGNP/360.0*2.0*PI;
	d = deltaGNP/360.0*2.0*PI;
	a = alphaGNP/360.0*2.0*PI;
	
	T[0][0] = -cos(t)*sin(d)*cos(a)-sin(t)*sin(a);
	T[0][1] = -cos(t)*sin(d)*sin(a)+sin(t)*cos(a);
	T[0][2] = cos(t)*cos(d);
	
	T[1][0] = -sin(t)*sin(d)*cos(a)+cos(t)*sin(a);
	T[1][1] = -sin(t)*sin(d)*sin(a)-cos(t)*cos(a);
	T[1][2] = sin(t)*cos(d);
	
	T[2][0] = cos(d)*cos(a);
	T[2][1] = cos(d)*sin(a);
	T[2][2] = sin(d);
	
	//invert matrix T in the most general way
	detT = T[0][0]*T[1][1]*T[2][2] + T[1][0]*T[2][1]*T[0][2] + T[2][0]*T[0][1]*T[1][2] - T[0][0]*T[1][2]*T[2][1] - T[1][0]*T[2][2]*T[0][1] - T[2][0]*T[0][2]*T[1][1];
	
	TI[0][0] = (T[1][1]*T[2][2]-T[1][2]*T[2][1])/detT;
	TI[1][0] = (T[1][2]*T[2][0]-T[1][0]*T[2][2])/detT;
	TI[2][0] = (T[1][0]*T[2][1]-T[1][1]*T[2][0])/detT;
	
	TI[0][1] = (T[0][2]*T[2][1]-T[0][1]*T[2][2])/detT;
	TI[1][1] = (T[0][0]*T[2][2]-T[2][0]*T[0][2])/detT;
	TI[2][1] = (T[0][1]*T[2][0]-T[0][0]*T[2][1])/detT;
	
	TI[0][2] = (T[0][1]*T[1][2]-T[0][2]*T[1][1])/detT;
	TI[1][2] = (T[0][2]*T[1][0]-T[0][0]*T[1][2])/detT;
	TI[2][2] = (T[0][0]*T[1][1]-T[0][1]*T[1][0])/detT;
	
	//convert to kpc
	x = xtemp[0]/1000.0;
	y = xtemp[1]/1000.0;
	z = xtemp[2]/1000.0;

	dsun = *dsuntemp/1000.0;
	
	vx = vtemp[0];
	vy = vtemp[1];
	vz = vtemp[2];

	vr = *vrtemp;
	vrsun = *vrsuntemp;
	
	//convert to radians
	DECrad = *DECtemp/360.0*2.0*PI;
	RArad = *RAtemp/360.0*2.0*PI;
	PArad = *PAtemp/360.0*2.0*PI;
	
	
	//get the galactic coordinates first
	if (coordtype == 1) {
		if (radio) printf("\nConverting equatorial to galactic coordinates using the transformation matrix:\n"); 
		if (radio) printf("%f\t%f\t%f\n",T[0][0],T[0][1],T[0][2]);
		if (radio) printf("%f\t%f\t%f\n",T[1][0],T[1][1],T[1][2]);
		if (radio) printf("%f\t%f\t%f\n",T[2][0],T[2][1],T[2][2]);
		
		brad = asin(T[2][0]*cos(DECrad)*cos(RArad) + T[2][1]*cos(DECrad)*sin(RArad) + T[2][2]*sin(DECrad));
		if (asin((T[1][0]*cos(DECrad)*cos(RArad) + T[1][1]*cos(DECrad)*sin(RArad) + T[1][2]*sin(DECrad))/cos(brad))>=0.0) {
			lrad = acos((T[0][0]*cos(DECrad)*cos(RArad) + T[0][1]*cos(DECrad)*sin(RArad) + T[0][2]*sin(DECrad))/cos(brad)); 
		} else {
			lrad = 2.0*PI-acos((T[0][0]*cos(DECrad)*cos(RArad) + T[0][1]*cos(DECrad)*sin(RArad) + T[0][2]*sin(DECrad))/cos(brad)); 			
		}
		lcosbrad = lrad*cos(brad);
	} else if (coordtype == 2) {
		brad = *btemp/360.0*2.0*PI;
		if (*ltemp) {
			lrad = *ltemp/360.0*2.0*PI;
			lcosbrad = lrad*cos(brad);
		} else if (*lcosbtemp) {
			lcosbrad = *lcosbtemp/360.0*2.0*PI;
			lrad = lcosbrad/cos(brad);
		}
	} else if (coordtype == 3) {
		if (y >= 0.0)
			lrad = acos((x-xsun)/sqrt(pow(x-xsun,2)+y*y));
		else
			lrad = 2.0*Pi-acos((x-xsun)/sqrt(pow(x-xsun,2)+y*y));
		brad =  atan(z/sqrt(pow(x-xsun,2)+y*y));
		lcosbrad = lrad*cos(brad);
	}
	
	
	//get 3d position of cluster [kpc] from galactic coordinates
	if (coordtype < 3) {
		z = sin(brad)*dsun;
		dxy = sqrt(dsun*dsun-z*z);
		x = cos(lrad)*dxy + xsun;
		y = sin(lrad)*dxy;
	} else if (coordtype == 3) {
		dsun = sqrt(pow(x-xsun,2)+y*y+z*z);		
	}
	
	
	//finally get the equatorial coordinates from galactic coordinates
	if (coordtype > 1) {
		if (radio) printf("\nConverting galactic to equatorial coordinates using the transformation matrix:\n"); 
		if (radio) printf("%f\t%f\t%f\n",TI[0][0],TI[0][1],TI[0][2]);
		if (radio) printf("%f\t%f\t%f\n",TI[1][0],TI[1][1],TI[1][2]);
		if (radio) printf("%f\t%f\t%f\n",TI[2][0],TI[2][1],TI[2][2]);

		if (radio) {
			//unit matrix B = T * TI
			B[0][0] = T[0][0]*TI[0][0] + T[0][1]*TI[1][0] + T[0][2]*TI[2][0];
			B[0][1] = T[0][0]*TI[0][1] + T[0][1]*TI[1][1] + T[0][2]*TI[2][1];
			B[0][2] = T[0][0]*TI[0][2] + T[0][1]*TI[1][2] + T[0][2]*TI[2][2];
	
			B[1][0] = T[1][0]*TI[0][0] + T[1][1]*TI[1][0] + T[1][2]*TI[2][0];
			B[1][1] = T[1][0]*TI[0][1] + T[1][1]*TI[1][1] + T[1][2]*TI[2][1];
			B[1][2] = T[1][0]*TI[0][2] + T[1][1]*TI[1][2] + T[1][2]*TI[2][2];
	
			B[2][0] = T[2][0]*TI[0][0] + T[2][1]*TI[1][0] + T[2][2]*TI[2][0];
			B[2][1] = T[2][0]*TI[0][1] + T[2][1]*TI[1][1] + T[2][2]*TI[2][1];
			B[2][2] = T[2][0]*TI[0][2] + T[2][1]*TI[1][2] + T[2][2]*TI[2][2];
	
			printf("\nCalculating T*T^{-1} = 1 for consistency check:\n");
			printf("%f\t%f\t%f\n",B[0][0],B[0][1],B[0][2]);
			printf("%f\t%f\t%f\n",B[1][0],B[1][1],B[1][2]);
			printf("%f\t%f\t%f\n",B[2][0],B[2][1],B[2][2]);
		}
		
		DECrad = asin(TI[2][0]*cos(brad)*cos(lrad)+TI[2][1]*cos(brad)*sin(lrad)+TI[2][2]*sin(brad));
		if (asin((TI[1][0]*cos(brad)*cos(lrad) + TI[1][1]*cos(brad)*sin(lrad) + TI[1][2]*sin(brad))/cos(DECrad))>=0.0) {
			RArad = acos((TI[0][0]*cos(brad)*cos(lrad) + TI[0][1]*cos(brad)*sin(lrad) + TI[0][2]*sin(brad))/cos(DECrad)); 
		} else {
			RArad = 2.0*PI-acos((TI[0][0]*cos(brad)*cos(lrad) + TI[0][1]*cos(brad)*sin(lrad) + TI[0][2]*sin(brad))/cos(DECrad)); 			
		}
	}
	
	

	//get tangential velocity in [km/s] from different sets of velocity-measurement types
	
	//get coordinates of equatorial north pole on great circle
	bENP = asin(T[2][0]*cos(DECENP)*cos(RAENP) + T[2][1]*cos(DECENP)*sin(RAENP) + T[2][2]*sin(DECENP));
	if (asin((T[1][0]*cos(DECENP)*cos(RAENP) + T[1][1]*cos(DECENP)*sin(RAENP) + T[1][2]*sin(DECENP))/cos(bENP))>=0.0) {
		lENP = acos((T[0][0]*cos(DECENP)*cos(RAENP) + T[0][1]*cos(DECENP)*sin(RAENP) + T[0][2]*sin(DECENP))/cos(bENP)); 
	} else {
		lENP = 2.0*PI-acos((T[0][0]*cos(DECENP)*cos(RAENP) + T[0][1]*cos(DECENP)*sin(RAENP) + T[0][2]*sin(DECENP))/cos(bENP)); 			
	}
	if (radio) printf("\nCoordinates of equatorial north pole:\n");
	if (radio) printf("bENP = %f\tlENP = %f\n", bENP, lENP);
	zENP = sin(bENP)*dsun;
	dxyENP = sqrt(dsun*dsun-zENP*zENP);
	xENP = cos(lENP)*dxyENP + xsun;
	yENP = sin(lENP)*dxyENP;
	if (radio) printf("xENP = %f\tyENP = %f\tzENP = %f\n", xENP, yENP, zENP);
	
	
	if (vcoordtype == 1) {

		//get radial velocity in 3d coordinates [km/s]
		vrx = (x - xsun)/dsun*vrsun;
		vry = y/dsun*vrsun;
		vrz = z/dsun*vrsun;
		if (radio) printf("\nHeliocentric radial velocity in cartesian coordinates:\n");
		if (radio) printf("vrx = %.3f\tvry = %.3f\tvrz = %.3f\tvr = %.3f [km/s] (heliocentric)\n",vrx,vry,vrz,sqrt(vrx*vrx+vry*vry+vrz*vrz));		
		
		//convert to km/s
		mu = *mutemp*dsun*4.74057;

		//compute proper motion components
		mu_alphacosdelta = mu*sin(PArad);
		mu_delta = mu*cos(PArad);
		
		A[0][0] = cos(RArad)*cos(DECrad);
		A[0][1] = -sin(RArad);
		A[0][2] = -cos(RArad)*sin(DECrad);
		
		A[1][0] = sin(RArad)*cos(DECrad);
		A[1][1] = cos(RArad);
		A[1][2] = -sin(RArad)*sin(DECrad);
		
		A[2][0] = sin(DECrad);
		A[2][1] = 0.0;
		A[2][2] = cos(DECrad);
		
		//printf("%f\t%f\t%f\n",A[0][0],A[0][1],A[0][2]);
		//printf("%f\t%f\t%f\n",A[1][0],A[1][1],A[1][2]);
		//printf("%f\t%f\t%f\n",A[2][0],A[2][1],A[2][2]);
		
		//B = T * A
		B[0][0] = T[0][0]*A[0][0] + T[0][1]*A[1][0] + T[0][2]*A[2][0];
		B[0][1] = T[0][0]*A[0][1] + T[0][1]*A[1][1] + T[0][2]*A[2][1];
		B[0][2] = T[0][0]*A[0][2] + T[0][1]*A[1][2] + T[0][2]*A[2][2];
		
		B[1][0] = T[1][0]*A[0][0] + T[1][1]*A[1][0] + T[1][2]*A[2][0];
		B[1][1] = T[1][0]*A[0][1] + T[1][1]*A[1][1] + T[1][2]*A[2][1];
		B[1][2] = T[1][0]*A[0][2] + T[1][1]*A[1][2] + T[1][2]*A[2][2];
		
		B[2][0] = T[2][0]*A[0][0] + T[2][1]*A[1][0] + T[2][2]*A[2][0];
		B[2][1] = T[2][0]*A[0][1] + T[2][1]*A[1][1] + T[2][2]*A[2][1];
		B[2][2] = T[2][0]*A[0][2] + T[2][1]*A[1][2] + T[2][2]*A[2][2];
		
		//printf("%f\t%f\t%f\n",B[0][0],B[0][1],B[0][2]);
		//printf("%f\t%f\t%f\n",B[1][0],B[1][1],B[1][2]);
		//printf("%f\t%f\t%f\n",B[2][0],B[2][1],B[2][2]);
		
		vx = vrsun*B[0][0] + mu_alphacosdelta*B[0][1] + mu_delta*B[0][2] +vxsun;
		vy = vrsun*B[1][0] + mu_alphacosdelta*B[1][1] + mu_delta*B[1][2] +vysun+vLSRtemp;
		vz = vrsun*B[2][0] + mu_alphacosdelta*B[2][1] + mu_delta*B[2][2] +vzsun;
		
		
		if (radio) printf("\nCartesian velocity:\n");
		if (radio) printf("vx = %.3f\tvy = %.3f\tvz = %.3f\tv = %.3f [km/s]\n",vx,vy,vz, sqrt(vx*vx+vy*vy+vz*vz));
		if (radio) printf("vx = %.3f\tvy = %.3f\tvz = %.3f\tv = %.3f [km/s] (heliocentric)\n",vx-vxsun,vy-vysun-vLSRtemp,vz-vzsun, sqrt(pow(vx-vxsun,2)+pow(vy-vysun-vLSRtemp,2)+pow(vz-vzsun,2)));
		
        
        //get proper motion in galactic coordinates (Poleski 2013)
        C1 = sin(d)*cos(DECrad)-cos(d)*sin(DECrad)*cos(RArad-a);
        C2 = cos(d)*sin(RArad-a);
        mu_lcosb = (C1*mu_alphacosdelta+C2*mu_delta)/cos(brad);
        mu_b = (-C2*mu_alphacosdelta+C1*mu_delta)/cos(brad);

        
	} else if (vcoordtype == 2) {

		//get radial velocity in 3d coordinates [km/s]
		vrx = (x - xsun)/dsun*vrsun;
		vry = y/dsun*vrsun;
		vrz = z/dsun*vrsun;
		if (radio) printf("\nHeliocentric radial velocity in cartesian coordinates:\n");
		if (radio) printf("vrx = %.3f\tvry = %.3f\tvrz = %.3f\tvr = %.3f [km/s] (heliocentric)\n",vrx,vry,vrz,sqrt(vrx*vrx+vry*vry+vrz*vrz));		
		
		if (*mu_alphatemp) *mu_alphacosdeltatemp = *mu_alphatemp*cos(DECrad);
		else if (*mu_alphacosdeltatemp) *mu_alphatemp = *mu_alphacosdeltatemp/cos(DECrad);

		//convert to km/s
		mu_alphacosdelta = *mu_alphacosdeltatemp*dsun*4.74057;
		mu_delta = *mu_deltatemp*dsun*4.74057;
		mu = sqrt(mu_alphacosdelta*mu_alphacosdelta+mu_delta*mu_delta);
		
		A[0][0] = cos(RArad)*cos(DECrad);
		A[0][1] = -sin(RArad);
		A[0][2] = -cos(RArad)*sin(DECrad);
		
		A[1][0] = sin(RArad)*cos(DECrad);
		A[1][1] = cos(RArad);
		A[1][2] = -sin(RArad)*sin(DECrad);
		
		A[2][0] = sin(DECrad);
		A[2][1] = 0.0;
		A[2][2] = cos(DECrad);
		
		//printf("%f\t%f\t%f\n",A[0][0],A[0][1],A[0][2]);
		//printf("%f\t%f\t%f\n",A[1][0],A[1][1],A[1][2]);
		//printf("%f\t%f\t%f\n",A[2][0],A[2][1],A[2][2]);
		
		//B = T * A
		B[0][0] = T[0][0]*A[0][0] + T[0][1]*A[1][0] + T[0][2]*A[2][0];
		B[0][1] = T[0][0]*A[0][1] + T[0][1]*A[1][1] + T[0][2]*A[2][1];
		B[0][2] = T[0][0]*A[0][2] + T[0][1]*A[1][2] + T[0][2]*A[2][2];
		
		B[1][0] = T[1][0]*A[0][0] + T[1][1]*A[1][0] + T[1][2]*A[2][0];
		B[1][1] = T[1][0]*A[0][1] + T[1][1]*A[1][1] + T[1][2]*A[2][1];
		B[1][2] = T[1][0]*A[0][2] + T[1][1]*A[1][2] + T[1][2]*A[2][2];
		
		B[2][0] = T[2][0]*A[0][0] + T[2][1]*A[1][0] + T[2][2]*A[2][0];
		B[2][1] = T[2][0]*A[0][1] + T[2][1]*A[1][1] + T[2][2]*A[2][1];
		B[2][2] = T[2][0]*A[0][2] + T[2][1]*A[1][2] + T[2][2]*A[2][2];
		
		//printf("%f\t%f\t%f\n",B[0][0],B[0][1],B[0][2]);
		//printf("%f\t%f\t%f\n",B[1][0],B[1][1],B[1][2]);
		//printf("%f\t%f\t%f\n",B[2][0],B[2][1],B[2][2]);
		
		vx = vrsun*B[0][0] + mu_alphacosdelta*B[0][1] + mu_delta*B[0][2] +vxsun;
		vy = vrsun*B[1][0] + mu_alphacosdelta*B[1][1] + mu_delta*B[1][2] +vysun+vLSRtemp;
		vz = vrsun*B[2][0] + mu_alphacosdelta*B[2][1] + mu_delta*B[2][2] +vzsun;
				
		if (radio) printf("\nCartesian velocity:\n");
		if (radio) printf("vx = %.3f\tvy = %.3f\tvz = %.3f\tv = %.3f [km/s]\n",vx,vy,vz, sqrt(vx*vx+vy*vy+vz*vz));
		if (radio) printf("vx = %.3f\tvy = %.3f\tvz = %.3f\tv = %.3f [km/s] (heliocentric)\n",vx-vxsun,vy-vysun-vLSRtemp,vz-vzsun, sqrt(pow(vx-vxsun,2)+pow(vy-vysun-vLSRtemp,2)+pow(vz-vzsun,2)));

		//get position angle of proper motion
		
		//heliocentric transverse velocity
		vtx = vx-vxsun-vrx;
		vty = vy-vysun-vLSRtemp-vry;
		vtz = vz-vzsun-vrz;
		if (radio) printf("\nTransverse velocity:\n");
		if (radio) printf("vtx = %f\tvty = %f\tvtz = %f\tvt = %f [km/s] (heliocentric)\n", vtx, vty, vtz, sqrt(vtx*vtx+vty*vty+vtz*vtz));
		
		//get tangential vector pointing to ENP
		FAK = -((xENP-xsun)*(x-xsun)+yENP*y+zENP*z)/(pow(x-xsun,2)+y*y+z*z);
		xdelta = FAK*(x-xsun)+(xENP-xsun);
		ydelta = FAK*y+yENP;
		zdelta = FAK*z+zENP;
		
		//determine distance (pos or neg) of Xobject + Vt from plane connecting ENP, Xobject and observer for position angle
		nx = y*zENP-z*yENP;
		ny = z*(xENP-xsun)-(x-xsun)*zENP;
		nz = (x-xsun)*yENP-y*(xENP-xsun);
		dvt = nx*(x+vtx)+ny*(y+vty)+nz*(z+vtz)-nx*xsun;

		//get position angle of proper motion with respect to tangential vector pointing to ENP
		if (dvt <= 0) 
			PArad = acos((xdelta*vtx+ydelta*vty+zdelta*vtz)/(sqrt(vtx*vtx+vty*vty+vtz*vtz)*sqrt(xdelta*xdelta+ydelta*ydelta+zdelta*zdelta)));
		else 
			PArad = 2.0*PI-acos((xdelta*vtx+ydelta*vty+zdelta*vtz)/(sqrt(vtx*vtx+vty*vty+vtz*vtz)*sqrt(xdelta*xdelta+ydelta*ydelta+zdelta*zdelta)));
		
		if (radio) printf("\nProper motion and position angle:\n");
		if (radio) printf("mu = %f\tPA = %f\n", mu, PArad);

		
        //get proper motion in galactic coordinates (Poleski 2013)
        C1 = sin(d)*cos(DECrad)-cos(d)*sin(DECrad)*cos(RArad-a);
        C2 = cos(d)*sin(RArad-a);
        mu_lcosb = (C1*mu_alphacosdelta+C2*mu_delta)/cos(brad);
        mu_b = (-C2*mu_alphacosdelta+C1*mu_delta)/cos(brad);

        
	} else if (vcoordtype == 3) {

		if (radio) printf("\nCartesian velocity:\n");
		if (radio) printf("vx = %.3f\tvy = %.3f\tvz = %.3f\tv = %.3f [km/s]\n",vx,vy,vz, sqrt(vx*vx+vy*vy+vz*vz));
		if (radio) printf("vx = %.3f\tvy = %.3f\tvz = %.3f\tv = %.3f [km/s] (heliocentric)\n",vx-vxsun,vy-vysun-vLSRtemp,vz-vzsun, sqrt(pow(vx-vxsun,2)+pow(vy-vysun-vLSRtemp,2)+pow(vz-vzsun,2)));		
		
		//heliocentric radial velocity
		vrsun = ((vx-vxsun)*(x-xsun)+(vy-vysun-vLSRtemp)*y+(vz-vzsun)*z)/sqrt(pow(x-xsun,2)+y*y+z*z);

		//get radial velocity in 3d coordinates [km/s]
		vrx = (x - xsun)/dsun*vrsun;
		vry = y/dsun*vrsun;
		vrz = z/dsun*vrsun;
		if (radio) printf("\nHeliocentric radial velocity in cartesian coordinates:\n");
		if (radio) printf("vrx = %.3f\tvry = %.3f\tvrz = %.3f\tvr = %.3f [km/s] (heliocentric)\n",vrx,vry,vrz,sqrt(vrx*vrx+vry*vry+vrz*vrz));				
		
		//get position angle of proper motion
		
		//heliocentric transverse velocity
		vtx = vx-vxsun-vrx;
		vty = vy-vysun-vLSRtemp-vry;
		vtz = vz-vzsun-vrz;
		if (radio) printf("\nTransverse velocity:\n");
		if (radio) printf("vtx = %f\tvty = %f\tvtz = %f\tvt = %f [km/s] (heliocentric)\n", vtx, vty, vtz, sqrt(vtx*vtx+vty*vty+vtz*vtz));
		
		//get tangential vector pointing to ENP
		FAK = -((xENP-xsun)*(x-xsun)+yENP*y+zENP*z)/(pow(x-xsun,2)+y*y+z*z);
		xdelta = FAK*(x-xsun)+(xENP-xsun);
		ydelta = FAK*y+yENP;
		zdelta = FAK*z+zENP;
		
		//determine distance (pos or neg) of Xobject + Vt from plane connecting ENP, Xobject and observer for position angle
		nx = y*zENP-z*yENP;
		ny = z*(xENP-xsun)-(x-xsun)*zENP;
		nz = (x-xsun)*yENP-y*(xENP-xsun);
		dvt = nx*(x+vtx)+ny*(y+vty)+nz*(z+vtz)-nx*xsun;
		
		//get position angle of proper motion with respect to tangential vector pointing to ENP
		if (dvt <= 0) 
			PArad = acos((xdelta*vtx+ydelta*vty+zdelta*vtz)/(sqrt(vtx*vtx+vty*vty+vtz*vtz)*sqrt(xdelta*xdelta+ydelta*ydelta+zdelta*zdelta)));
		else 
			PArad = 2.0*PI-acos((xdelta*vtx+ydelta*vty+zdelta*vtz)/(sqrt(vtx*vtx+vty*vty+vtz*vtz)*sqrt(xdelta*xdelta+ydelta*ydelta+zdelta*zdelta)));
		
		if (radio) printf("\nProper motion and position angle:\n");
		if (radio) printf("mu = %f\tPA = %f\n", mu, PArad);
		
		mu = sqrt(vtx*vtx+vty*vty+vtz*vtz);
		mu_delta = mu*cos(PArad);
		mu_alphacosdelta = mu*sin(PArad);
		
        
        //get proper motion in galactic coordinates (Poleski 2013)
        C1 = sin(d)*cos(DECrad)-cos(d)*sin(DECrad)*cos(RArad-a);
        C2 = cos(d)*sin(RArad-a);
        mu_lcosb = (C1*mu_alphacosdelta+C2*mu_delta)/cos(brad);
        mu_b = (-C2*mu_alphacosdelta+C1*mu_delta)/cos(brad);

	}
	
	
	if (radio) printf("\nProper motion:\n");
	if (radio) printf("mu_alphacosdelta  = %f\tmu_delta = %f\tmu = %f [km/s]\t PA = %f\n", mu_alphacosdelta, mu_delta, mu, PArad);
		
	vr = (vx*(x-xsun)+vy*y+vz*z)/sqrt(pow(x-xsun,2)+y*y+z*z);
	if (radio) printf("\nRadial velocity:\n");
	if (radio) printf("vr = %.3f\tvr = %.3f (heliocentric) [km/s]\n", vr, vrsun);

	//consistency check with formula for GSR radial velocity from script of Steven Majewski
	if (radio) vrLSR = vrsun + (vxsun*cos(brad)*cos(lrad)+vysun*cos(brad)*sin(lrad)+vzsun*sin(brad));
	if (radio) vrGSR = vrLSR + vLSRtemp*cos(brad)*sin(lrad);
	if (radio) printf("\nConsistency check with formula for Galactic standard of rest (GSR) radial velocity (should be equal to vr):\n");
	if (radio) printf("vr_LSR = %f\tvr_GSR = %.3f [km/s]\n", vrLSR, vrGSR);
	
	
	
	//convert back to input units and write to output
	*xtemp = 1000.0*x;
	*(xtemp+1) = 1000.0*y;
	*(xtemp+2) = 1000.0*z;
	
	*dsuntemp = 1000.0*dsun;
	
	*vtemp = vx;
	*(vtemp+1) = vy;
	*(vtemp+2) = vz;

	*vrsuntemp = vrsun;
	*vrtemp = vr;

	*DECtemp = DECrad*180.0/PI;
	*RAtemp = RArad*180.0/PI;
	
	*btemp = brad*180.0/PI;
	*ltemp = lrad*180.0/PI;
	*lcosbtemp = *ltemp*cos(brad);

	*mutemp = mu/(dsun*4.74057);
	*PAtemp = PArad*180.0/PI;
	*mu_deltatemp = mu_delta/(dsun*4.74057);
	*mu_alphacosdeltatemp = mu_alphacosdelta/(dsun*4.74057);
	*mu_alphatemp = *mu_alphacosdeltatemp/cos(DECrad);
    
    *mu_btemp = mu_b/(dsun*4.74057);
    *mu_lcosbtemp = mu_lcosb/(dsun*4.74057);
    *mu_ltemp = *mu_lcosbtemp/cos(brad);

	
}
