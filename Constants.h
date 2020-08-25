// All of the defined parameters.  (Most are not used.)

#define SECSYR 3.15581497632e7  // Number of seconds in a sidereal year 
#define SECSMTH 2.629743834e6   // Number of seconds in a month     
#define SECSDAY 86400.0          // Number of seconds in a day            
#define CLIGHT 2.99792458e8     // Speed of light in m/s      
#define G 6.67259e-11           // Newton's constant  
#define TSUN  4.92549232189886339689643862e-6  // mass to seconds conversion   
#define AU 1.4959787e11         // Astronomical Unit in meters
#define EULER 0.577215664901532860606512  // Euler's constant           
#define MSOLAR 1.98892e30         // Mass of the sun in kg     
#define ln2 0.693147180559945   // ln 2
#define ln32 0.405465108108164  // ln(3/2)
#define SQ3 1.732050807568877  // sqrt(3)
#define PI 3.1415926535897931159979634685442      // Pi                        
#define TPI 6.2831853071795862319959269370884     // 2 Pi                      
#define PI2 9.869604401089357992304940125905      // Pi^2
#define PIon2 1.5707963267948966                   // PI/2
#define PIon3 1.0471975511965977                  //PI/3
#define IPI 0.31830988618379069121644420192752    // 1 / Pi                    
#define RPI 1.7724538509055158819194275565678     // sqrt(Pi)                  
#define IRPI 0.56418958354775627928034964497783   // 1 / sqrt(Pi)              
#define RTPI 2.5066282746310002416123552393401    // sqrt(2 Pi)                
#define IRTPI 0.39894228040143270286321808271168  // 1 / sqrt(2 Pi)
#define IRTWO 0.707106781186547549 // 1/sqrt(2)
#define PARSEC 3.08568025e16    // Parsec in meters          
#define GPC 3.08568025e25       // Gigaparsec in meters     
#define RAU 499.0047815         // 1 AU in seconds           
#define HUBBLE 2.30095130563e-18 // Hubble constant 71 km/s/Mpc
#define LMF 3.16887646154e-8    // LISA modulation frequency              


// Cash-Karp Runge-Kutta parameters:

#define A2 (1./5.)
#define A3 (3./10.)
#define A4 (3./5.)
#define A5 1.
#define A6 (7./8.)

#define B21 (1./5.)
#define B31 (3./40.)
#define B32 (9./40.)
#define B41 (3./10.)
#define B42 (-9./10.)
#define B43 (6./5.)
#define B51 (-11./54.)
#define B52 (5./2.)
#define B53 (-70./27.)
#define B54 (35./27.)
#define B61 (1631./55296.)
#define B62 (175./512.)
#define B63 (575./13824.)
#define B64 (44275./110592.)
#define B65 (253./4096.)

#define C1 (37./378.)
#define C2 0.
#define C3 (250./621.)
#define C4 (125./594.)
#define C5 0.
#define C6 (512./1771.)

#define D1 (2825./27648.)
#define D2 0.
#define D3 (18575./48384.)
#define D4 (13525./55296.)
#define D5 (277./14336.)
#define D6 (1./4.)

#define RK_H 100000
#define RK_EPS 0.00000001

// Are we on the last step of an update?

#define ENDYES 1
#define ENDNO 2

// Length of RK storage vectors:
#define VECLENGTH 10000
