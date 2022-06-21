#define  PHYSICS                        HD
#define  DIMENSIONS                     3
#define  GEOMETRY                       CARTESIAN
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 PARABOLIC
#define  TIME_STEPPING                  CHARACTERISTIC_TRACING
#define  NTRACER                        2
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            10

/* -- physics dependent declarations -- */

#define  DUST_FLUID                     NO
#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO

/* -- user-defined parameters (labels) -- */

#define  PAR_POWER                      0
#define  PAR_MACH                       1
#define  PAR_SPEED                      2
#define  PAR_RADIUS                     3
#define  PAR_RHOW                       4
#define  PAR_RHOH                       5
#define  PAR_THOT                       6
#define  PAR_TCRIT                      7
#define  PAR_RADIUS_W                   8
#define  PAR_SIGMA                      9

/* [Beg] user-defined constants (do not change this line) */

#define  UNIT_DENSITY                   (0.6036*CONST_amu)
#define  UNIT_VELOCITY                  CONST_c
#define  UNIT_LENGTH                    (1000*CONST_pc)
#define  START_MODE_AGN                 0
#define  START_MODE_CLOUD               1
#define  START_MODE_HALO                2
#define  START_MODE                     START_MODE_HALO




/* [End] user-defined constants (do not change this line) */
