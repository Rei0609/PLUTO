#define  PHYSICS                        RHD
#define  DIMENSIONS                     2
#define  GEOMETRY                       SPHERICAL
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  NTRACER                        0
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            4

/* -- physics dependent declarations -- */

#define  EOS                            TAUB
#define  ENTROPY_SWITCH                 NO
#define  RADIATION                      NO

/* -- user-defined parameters (labels) -- */

#define  BETA                           0
#define  RHO_IN                         1
#define  RHO_OUT                        2
#define  PRESS_IN                       3

/* [Beg] user-defined constants (do not change this line) */

#define  LIMITER                        OSPRE_LIM
#define  CHOMBO_LOGR                    YES
#define  RECONSTRUCT_4VEL               NO

/* [End] user-defined constants (do not change this line) */
