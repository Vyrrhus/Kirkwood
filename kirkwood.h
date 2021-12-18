#ifndef KIRKWOOD_H
    #define KIRKWOOD_H
    
    /***** SIMULATION PARAMETERS */
    #define ASTEROID_NUMBER 10000
    #define JUPITER_REV     10000

    #define STEP_TIME       0.1
    #define STEP_WRITING    1.
    #define SAVE_BOOL       0       // True (1) or False (0)
    #define SAVE_HIST       0       // Idem

    /***** ASTEROID SEMI MAJOR AXIS GENERATION */
    #define MIN_AST_GENERATION  1.6
    #define MAX_AST_GENERATION  5.2
    #define RESET_ASTEROID_INPUT 1     // True (1) or False (0)

    /***** DELETE CONDITIONS */
    // Asteroids with Kepler elements greater than these won't be computed anymore
    #define SMA_MAX_VALUE   15
    #define ECC_MAX_VALUE   1

    /***** ALIAS */
    #define AU_CONST        149597870700
    #define YEAR_CONST      86400 * 365.2524
    #define GM_CONVERSION   pow(YEAR_CONST,2) / pow(AU_CONST,3)
    #define JUPITER_YEAR    11.87
    #define STD_SUN         1.32712440018e+020 * GM_CONVERSION
    #define RAD_CONVERSION  M_PI/180

    #define LENGTH_SIMULATION   JUPITER_YEAR * JUPITER_REV
    // #define LENGTH_SIMULATION   1 * STEP_TIME


    /***** STRUCTURES */
    typedef struct {
        double x,y,z;
    } Vector;

    typedef struct {
        int stop_simu;
        Vector pos;
        Vector vel;
        Vector acc;
        double std;
        double semiMajorAxis;
        double eccentricity;
        double inclination;
        double longitudeAN;
        double argPeriapsis;
    } Body;

#endif