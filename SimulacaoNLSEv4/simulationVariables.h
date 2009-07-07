

#define NPOINTS_T 	(512)
#define DELTAT	     0.1e-15 //0.1fs

#define NPOINTS_R 	(512)
#define DELTAR		2e-6 //100 um

#define TIMESCALE 	"s"
#define DISTANCESCALE 	"m"
//dont forget to change Nature's constants if changing the scales

#define VERSION "NLEEv4"


///////////////////////////////////////////////////
//Prototypes
void 	Propagate(double step);
void 	RunSimulation();


///////////////////////////////////////////////////////////
// Preprocessor code to avoid definition in various .cpp's
#ifdef __MAINCPP__ 
#define EXTERN 
///////////////////////////////////////////////////
//  Nature's constants. Careful with the units.
EXTERN  double 			lightVelocity = 3e8;	//light velocity
EXTERN  double 			c = lightVelocity;	//light velocity
EXTERN  double 		    m_electron = 9.109e-31;			// electron mass
EXTERN  double 		    q_electron = 1.609e-19;			// electron charge
EXTERN  double 		    epslon0 = 8.854e-12;			// electric permittivity
#else
#define EXTERN extern
///////////////////////////////////////////////////
//  save variables as above but w/o initialization
EXTERN  double 			lightVelocity;
EXTERN  double 			c;
EXTERN  double 		    m_electron;
EXTERN  double 		    q_electron;
EXTERN  double 		    epslon0;
#endif

///////////////////////////////////////////////////
//  Simulation Variables to be loaded from Config
EXTERN	double		zDistance;	     //total distance to propagate impulse
EXTERN	double		zOutputStep;	 //z distance output step
EXTERN	double 		Pmax; 	         //Pulse peak power
EXTERN	double 		Energy; 	     //Pulse Energy
EXTERN	double 		pulseT;	         //Temporal Width of our initial impulse (E)
EXTERN	double 		curvature;	             //Curvature of the beam
EXTERN  double 		spot;		 //gaussian beam diameter
EXTERN	double 		GVD;		//Group velocity dispersion at the central wavelength
EXTERN	double 		absorption;	//Absorption of the medium
EXTERN	double		n2;		    //nonlinear refractive index
EXTERN	double 		lambdaZero;	//central wavelength of our pulse
EXTERN	double 		desiredError;	//desired Error for adaptive error estimation
EXTERN  double 		boundaryRatio;		//=sizeof(boundaryLayer)/MAXR
EXTERN  double 		boundarySlope;		//slope of the absorption coeficient exponential
EXTERN  double 		boundaryI0;			//intensity of the absorption coeficient exponential
EXTERN  double 		n;                  // refractive index
EXTERN  double 		sigma;              // cross section for inverse bremsstrahlung
EXTERN  double      Ui;                 // ionization potential
EXTERN  double 		sigmaK;             // cross section for multiphoton ionization
EXTERN  double 		K;                  // multiphoton ionization order
EXTERN  double 		rho_neutral;        // neutral species density
EXTERN  double 		tau_r;              // electron decay time in the plasma
EXTERN  double 		tau_c;              // electron collision time
EXTERN  double 		p;              // gas pressure ratio  (p/1 atm)

///////////////////////////////////////////////////
// Simulation Variables to be calculated
//EXTERN	int			iteraction;
EXTERN	double complex		initialE[NPOINTS_T*NPOINTS_R];   //initial envelopes
EXTERN	double complex		E[NPOINTS_T*NPOINTS_R];          //
EXTERN	double complex		fftE[NPOINTS_T*NPOINTS_R];       // envelopes' transform
EXTERN	double      		rho[NPOINTS_T*NPOINTS_R];        // electron density
 
EXTERN	double complex		triIn[NPOINTS_R]; //for tridiagonal calculation
EXTERN	double complex		triLeft[NPOINTS_R]; //for tridiagonal calculation
EXTERN	double complex		triDiag[NPOINTS_R]; //for tridiagonal calculation
EXTERN	double complex		triRight[NPOINTS_R]; //for tridiagonal calculation
EXTERN	double complex		triOut[NPOINTS_R];  //for tridiagonal calculation
EXTERN	double complex		h2E[NPOINTS_T*NPOINTS_R]; //for adaptive calculation
EXTERN	double 				absorptionCalc[NPOINTS_R]; //for boundary layer 

EXTERN	double 			omegaZero;	//central frequency of our pulse
EXTERN  double 			rhoC;               // critical electron density



///////////////////////////////////////////////////
//  Maintenance Variables
EXTERN	fftw_plan    	forward[NPOINTS_R],backward[NPOINTS_R];
EXTERN	int 			startTime;
EXTERN	string 			configFileName;
EXTERN  string  		simDir;
EXTERN	struct tm * 	timeinfo;


