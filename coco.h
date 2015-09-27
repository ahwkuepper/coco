#define Pi 3.14159265
#define PI 3.14159265
#define G  0.0043009211           //gravitational constant in [km^2/s^2/Msun*pc]
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )


void convert(double *x, double *v, double *dsun, double *vrsun, double *vr, double *l, double *b, double *lcosb, double *RA, double *DEC, double *mu_alpha, double *mu_alphacosdelta, double *mu_delta, double *mutemp, double *PAtemp, double *mu_ltemp, double *mu_btemp, double *mu_lcosbtemp, int coordtype, int vcoordtype, double vLSRtemp);


//Galactic North Pole parameters
//double const alphaGNP = 192.25; //Galactic north pole in B1950 coordinates
//double const deltaGNP = 27.4;
//double const PAGNP = 123;       //Position angle with respect to equatorial pole
double const alphaGNP = 192.859508; //Galactic north pole in J2000 coordinates
double const deltaGNP = 27.128336;
double const PAGNP = 122.932;       //Position angle with respect to equatorial pole

//solar parameters
double const rgalsun = 8330.0; //solar Galactocentric radius [pc] (standard = 8330.0; Gillessen et al. 2009)
double const vLSR = 239.5;    //rotational velocity of local standard of rest

//solar reflex motion
//double const vxsun = 0.0;   //off
//double const vysun = 0.0;
//double const vzsun = 0.0;
//double const vxsun = 8.83; //±0.24 - solar motion with respect to the LSR from Coskunoglu et al. (2011) [km/s]
//double const vysun = 14.19; //±0.34 - 20453 RAVE stars
//double const const vzsun = 6.57; //±0.21
double const vxsun = 11.1;//+0.69/−0.75 - solar motion with respect to the LSR from Schoenrich, Binney & Dehnen (2010) [km/s]
double const vysun = 12.24;//+0.47−0.47
double const vzsun = 7.25;//+0.37−0.36
//double const vxsun = 10.0;    //solar motion with respect to the LSR from Dehnen & Binney (1998) [km/s]
//double const vysun = 5.3;
//double const vzsun = 7.2;
//double const vxsun = 10.4;    //solar motion with respect to the LSR from Johnson & Soderblom (1987) [km/s]
//double const vysun = 14.8;
//double const vzsun = 7.3;
//double const vxsun = 11.0;    //solar motion with respect to the LSR from Ratnatunga, Bahcall & Casertano (1989) [km/s]
//double const vysun = 14.0;
//double const vzsun = 7.5;


int const radio = 0; //say what you're doing (0 = no, 1 = yes)

