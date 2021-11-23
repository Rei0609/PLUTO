#define  PHYSICS                        HD
#define  DIMENSIONS                     1
#define  GEOMETRY                       SPHERICAL
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 PARABOLIC
#define  TIME_STEPPING                  RK3
#define  NTRACER                        1
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            8

/* -- physics dependent declarations -- */

#define  DUST_FLUID                     NO
#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 SELECTIVE
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO

/* -- user-defined parameters (labels) -- */

#define  PAR_POWER                      0
#define  PAR_MACH                       1
#define  PAR_SPEED                      2
#define  PAR_NHOT                       3
#define  PAR_THOT                       4
#define  PAR_TVIR                       5
#define  PAR_KAPPA                      6
#define  PAR_LAMBDA                     7

/* [Beg] user-defined constants (do not change this line) */

#define  UNIT_DENSITY                   (0.6036*CONST_amu)
#define  UNIT_VELOCITY                  CONST_c
#define  UNIT_LENGTH                    (1000.*CONST_pc)

/* [End] user-defined constants (do not change this line) */
