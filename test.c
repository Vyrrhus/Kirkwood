#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "engine.h"
#include "kirkwood.h"

/***** FUNCTIONS */
// Integrators
void leapfrog() {
    
};

// Physic model
void kepler_from_state(Body sun, Body *body) {
    /*  Compute kepler elements (a, e, w) from state vectors (position, velocity)
     */
    double a, e, w, energy;
    double sunStdParameter = (double) STD_SUN;

    // Relative state vectors to the Sun
    Vector relPos = {body->pos.x - sun.pos.x, 
                     body->pos.y - sun.pos.y};
    Vector relVel = {body->vel.x - sun.vel.x, 
                     body->vel.y - sun.vel.y};

    // Squared distance and velocity
    double r2 = relPos.x * relPos.x + relPos.y * relPos.y;
    double v2 = relVel.x * relVel.x + relVel.y * relVel.y;

    // Orbital momentum
    double h = relPos.x * relVel.y - relPos.y * relVel.x;

    // Eccentricity vector
    Vector ecc = {(relVel.y * h) / sunStdParameter - relPos.x / sqrt(r2),
                 -(relVel.x * h) / sunStdParameter - relPos.y / sqrt(r2)};
    
    // Mechanical specific energy, semi-major axis & eccentricity
    energy  = v2 / 2 - STD_SUN / sqrt(r2);
    a       = - STD_SUN / (2 * energy);
    e       = sqrt(ecc.x * ecc.x + ecc.y * ecc.y);



    // True longitude for near-circular orbit
    if (e < 1e-6) {
        w = acos(relPos.x / sqrt(r2));
        if (relPos.y > 0) {
            w = 2 * M_PI - w;
        }
    } else {
    // ... for non-circular orbit
        w = acos(ecc.x / e);
        if (ecc.y < 0) {
            w = 2 * M_PI - w;
        }
    }

    // If the body is the Sun, (a,e,w) = 0
    if (r2 != 0) {
        body->semiMajorAxis = a;
        body->eccentricity  = e;
        body->trueLongitude = w;
    }
};
void state_from_kepler(Body sun, Body *body) {
    /*  Compute state vectors (pos, vel) from Kepler elements (a, e, w).
     */
    double a, e, w, r, v;
    double cosw, sinw;

    a = body->semiMajorAxis;
    e = body->eccentricity;
    w = body->trueLongitude;

    cosw = cos(w);
    sinw = sin(w);

    // Distance & velocity relative to the Sun
    r = a * (1 - e);
    v = sqrt(STD_SUN * a * (1 - e * e)) / r;

    // State vectors
    body->pos.x =  r * cosw + sun.pos.x;
    body->pos.y = -r * sinw + sun.pos.y;
    body->vel.x =  v * sinw + sun.vel.x;
    body->vel.y =  v * cosw + sun.vel.y;
};










#ifndef ENGINE_H
    #define ENGINE_H
    
    #include "kirkwood.h"

    /***** PROTOTYPES */
    // Integrators
    void   leapfrog();
    void    yoshida();
    void yoshida8th();

    // Physic model
    void get_acceleration();
    void kepler_from_state(Body sun, Body *body);
    void state_from_kepler(Body sun, Body *body);

#endif



#include <stdio.h>
#include <stdlib.h>
#include "kirkwood.h"
#include "engine.h"



#ifndef KIRKWOOD_H
    #define KIRKWOOD_H
    
    /***** SIMULATION PARAMETERS */
    #define ASTEROID_NUMBER 10
    #define JUPITER_REV     100

    #define STEP_TIME       0.5
    #define STEP_WRITING    1.

    /***** ASTEROID GENERATION PARAMETERS */
    const double semiMajorAxisRange[2] = {1.9, 3.6};

    /***** ALIAS */
    #define AU_CONST        149597870700
    #define YEAR_CONST      86400 * 365.2524
    #define GM_CONVERSION   pow(YEAR_CONST,2) / pow(AU_CONST,3)
    #define JUPITER_YEAR    11.87
    #define STD_SUN         1.32712440018e+020 * GM_CONVERSION

    #define LENGTH_SIMULATION   JUPITER_YEAR * JUPITER_REV

    /***** STRUCTURES */
    typedef struct {
        double x,y;
    } Vector;

    typedef struct {
        Vector pos;
        Vector vel;
        Vector acc;
        double semiMajorAxis;
        double eccentricity;
        double trueLongitude;
    } Body;

#endif








