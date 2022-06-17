#define  PHYSICS                        HD
#define  DIMENSIONS                     1
#define  GEOMETRY                       CARTESIAN
#define  BODY_FORCE                     NO
#define  COOLING                        TABULATED
#define  RECONSTRUCTION                 PARABOLIC
#define  TIME_STEPPING                  RK3
#define  NTRACER                        1
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            5

/* -- physics dependent declarations -- */

#define  DUST_FLUID                     NO
#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 SELECTIVE
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO

/* -- user-defined parameters (labels) -- */

#define  PAR_RHO1                       0
#define  PAR_RHO2                       1
#define  PAR_TE2                        2
#define  PAR_P1P2                       3
#define  PAR_MACH                       4

/* [Beg] user-defined constants (do not change this line) */

#define  MU_NORM                        0.60364
#define  UNIT_DENSITY                   (CONST_amu*MU_NORM)
#define  UNIT_LENGTH                    (CONST_pc)
#define  UNIT_VELOCITY                  (1.e5)
#define  MU_CALC                        MU_ANALYTIC
#define  MU_CONST                       0
#define  MU_TABLE                       1
#define  MU_ANALYTIC                    2
#define  MU_FRACTIONS                   3
#define  MU_FNAME                       "mutable.dat"

/* [End] user-defined constants (do not change this line) */
