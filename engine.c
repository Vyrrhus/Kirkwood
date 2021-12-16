#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "engine.h"
#include "kirkwood.h"

/***** FUNCTIONS */
// Static functions
static Vector get_acceleration(Body body, int nbMajorBodies, Body listMajorBodies[]) {
    /*
     *  Compute the acceleration vector of a given body in a N-body system
     *  Main attractors are within listMajorBodies[]
     * 
     *  ----------------
     *      WARNING
     *      Collisions are not handled and might cause failure
     *      ie distance with attractors = 0
     */
    Vector acceleration = {0,0};
    Vector relativePosition;
    double r2;
    double a;

    // printf("   => FUNCTION\n");
    // printf("      Acc ini = {%g ; %g}\n", acceleration.x, acceleration.y);
    // printf("      Nb major bodies : {%d}\n", nbMajorBodies);
    // printf("      Position : {%g ; %g}\n", body.pos.x, body.pos.y);

    for (int i = 0 ; i < nbMajorBodies ; i++) {
        relativePosition.x = body.pos.x - listMajorBodies[i].pos.x;
        relativePosition.y = body.pos.y - listMajorBodies[i].pos.y;
        // printf("      =Body {%d} :  position {%g ; %g}\n",i, listMajorBodies[i].pos.x, listMajorBodies[i].pos.y);
        // printf("                   relative {%g ; %g}\n", relativePosition.x, relativePosition.y);
        // printf("                   std para {%g}\n", listMajorBodies[i].std);
        r2 = relativePosition.x * relativePosition.x + relativePosition.y * relativePosition.y;
        a  = listMajorBodies[i].std / (r2 * sqrt(r2));
        // printf("                   r2 = %g\n", r2);
        // printf("                   a  = %g\n", a);

        if (r2 != 0) {
            acceleration.x -= relativePosition.x * a;
            acceleration.y -= relativePosition.y * a;
                    // printf("                   accelera {%g ; %g}\n",acceleration.x, acceleration.y);

        }
    }

    return acceleration;
};

// Integrators
Body leapfrog(int nbMajorBodies, Body listMajorBodies[], Body listAsteroids[], double dT) {
    /*
     *  Leapfrog integrator
     */
    Vector acc = {0, 0};
    Body sun;

    // printf("\nLEAPFROG\n");

    // Position at t+dt
    // printf("ASTEROIDS\n");
    for (int i = 0 ; i < ASTEROID_NUMBER ; i++) {
        // printf("%d : %g, %g\n", i, listAsteroids[i].pos.x, listAsteroids[i].pos.y);
        // printf("---- %g, %g\n", listAsteroids[i].vel.x, listAsteroids[i].vel.y);
        // printf("---- %g, %g\n", listAsteroids[i].acc.x, listAsteroids[i].acc.y);

        listAsteroids[i].pos.x += listAsteroids[i].vel.x * dT + 0.5 * listAsteroids[i].acc.x * dT * dT;
        listAsteroids[i].pos.y += listAsteroids[i].vel.y * dT + 0.5 * listAsteroids[i].acc.y * dT * dT;
        // printf("=new: %g, %g\n", listAsteroids[i].pos.x, listAsteroids[i].pos.y);
    }
    // printf("MAJOR BODIES\n");
    for (int i = 0 ; i < nbMajorBodies ; i++) {
        // printf("%d : %g, %g\n", i, listMajorBodies[i].pos.x, listMajorBodies[i].pos.y);
        // printf("---- %g, %g\n", listMajorBodies[i].vel.x, listMajorBodies[i].vel.y);
        // printf("---- %g, %g\n", listMajorBodies[i].acc.x, listMajorBodies[i].acc.y);

        listMajorBodies[i].pos.x += listMajorBodies[i].vel.x * dT + 0.5 * listMajorBodies[i].acc.x * dT * dT;
        listMajorBodies[i].pos.y += listMajorBodies[i].vel.y * dT + 0.5 * listMajorBodies[i].acc.y * dT * dT;
        // printf("=new: %g, %g\n", listMajorBodies[i].pos.x, listMajorBodies[i].pos.y);
    }

    // Velocity and acceleration at t+dt
    // printf("ASTEROIDS\n");
    for (int i = 0 ; i < ASTEROID_NUMBER ; i++) {
        // printf("== nb {%d}", i);
        acc = get_acceleration(listAsteroids[i], nbMajorBodies, listMajorBodies);
        // printf("%d : %g, %g\n", i, listAsteroids[i].vel.x, listAsteroids[i].vel.y);
        // printf("   acc: {%g ; %g}\n", acc.x, acc.y);
        listAsteroids[i].vel.x += 0.5 * (listAsteroids[i].acc.x + acc.x) * dT;
        listAsteroids[i].vel.y += 0.5 * (listAsteroids[i].acc.y + acc.y) * dT;
        listAsteroids[i].acc = acc;
        // printf("---- %g, %g\n", listAsteroids[i].vel.x, listAsteroids[i].vel.y);
    }
    // printf("MAJOR BODIES\n");
    for (int i = 0 ; i < nbMajorBodies ; i++) {
        // printf("== nb {%d}", i);
        acc = get_acceleration(listMajorBodies[i], nbMajorBodies, listMajorBodies);
        // printf("%d : %g, %g\n", i, listMajorBodies[i].vel.x, listMajorBodies[i].vel.y);
        // printf("   acc: {%g ; %g}\n", acc.x, acc.y);
        listMajorBodies[i].vel.x += 0.5 * (listMajorBodies[i].acc.x + acc.x) * dT;
        listMajorBodies[i].vel.y += 0.5 * (listMajorBodies[i].acc.y + acc.y) * dT;
        listMajorBodies[i].acc = acc;
        // printf("---- %g, %g\n", listMajorBodies[i].vel.x, listMajorBodies[i].vel.y);

        if (listMajorBodies[i].std == STD_SUN) {
            sun = listMajorBodies[i];
        }
    }

    return sun;
};
Body yoshida(int nbMajorBodies, Body listMajorBodies[], Body listAsteroids[], double dT) {
    /*
     *  4th order Yoshida integrator
     */
    Vector acc = {0, 0};

    double a = pow(2, (double) 1/3);
    double W0 = - a / (2 - a);
    double W1 =   1 / (2 - a);
    double c[4] = {W1/2, (W0+W1)/2, (W0+W1)/2, W1/2};
    double d[3] = {W1, W0, W1};

    // Yoshida steps
    for (int step = 0 ; step < 4 ; step++) {
        // Position
        for (int i = 0 ; i < ASTEROID_NUMBER ; i++) {
            if (listAsteroids[i].isTooFar) {
                continue;
            }
            listAsteroids[i].pos.x += c[step] * listAsteroids[i].vel.x * dT;
            listAsteroids[i].pos.y += c[step] * listAsteroids[i].vel.y * dT;
        }
        for (int i = 0 ; i < nbMajorBodies ; i++) {
            listMajorBodies[i].pos.x += c[step] * listMajorBodies[i].vel.x * dT;
            listMajorBodies[i].pos.y += c[step] * listMajorBodies[i].vel.y * dT;
        }

        // No need to compute velocity for the 4th step
        if (step == 3) {
            break;
        }

        // Velocity and acceleration
        for (int i = 0 ; i < ASTEROID_NUMBER ; i++) {
            if (listAsteroids[i].isTooFar) {
                continue;
            }
            acc = get_acceleration(listAsteroids[i], nbMajorBodies, listMajorBodies);
            listAsteroids[i].vel.x += d[step] * acc.x * dT;
            listAsteroids[i].vel.y += d[step] * acc.y * dT;
            listAsteroids[i].acc = acc;
        }
        for (int i = 0 ; i < nbMajorBodies ; i++) {
            acc = get_acceleration(listMajorBodies[i], nbMajorBodies, listMajorBodies);
            listMajorBodies[i].vel.x += d[step] * acc.x * dT;
            listMajorBodies[i].vel.y += d[step] * acc.y * dT;
            listMajorBodies[i].acc = acc;
        }
    }

    // Return Sun (ie central body)
    for (int i = 0 ; i < nbMajorBodies ; i++) {
        if (listMajorBodies[i].std == STD_SUN) {
            return listMajorBodies[i];
        }
    }
};

// Physic model
void kepler_from_state(Body sun, Body *body) {
    /*  Compute kepler elements (a, e, w) from state vectors (position, velocity)
     */
    double a, e, w, energy;
    double sunStdParameter = (double) STD_SUN;

    // printf("BODY STD : {%g}\n", body->std);
    // printf("     pos: {%g ; %g}\n", body->pos.x, body->pos.y);
    // printf("     vel: {%g ; %g}\n", body->vel.x, body->vel.y);


    // Relative state vectors to the Sun
    Vector relPos = {body->pos.x - sun.pos.x, 
                     body->pos.y - sun.pos.y};
    Vector relVel = {body->vel.x - sun.vel.x, 
                     body->vel.y - sun.vel.y};

    // Squared distance and velocity
    double r2 = relPos.x * relPos.x + relPos.y * relPos.y;
    double v2 = relVel.x * relVel.x + relVel.y * relVel.y;

    // printf("== r2 = %g\n", r2);
    // printf("SUN pos: {%g ; %g}\n", sun.pos.x, sun.pos.y);
    // printf("BODYpos: {%g ; %g}\n", body->pos.x, body->pos.y);
    // Sun
    if (r2 == 0) {
        return;
    }

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

    body->semiMajorAxis = a;
    body->eccentricity  = e;
    body->trueLongitude = w;

    double rp, ra;
    rp = a * (1 - e);
    ra = a * (1 + e);

    if (rp <= 1.666 || ra >= 4.951) {
            body->eccentricity = 1;
            body->isTooFar = 1;
    }

    // printf("       a: {%g}\n", body->semiMajorAxis);
    // printf("       e: {%g}\n", body->eccentricity);
    // printf("       w: {%g}\n", body->trueLongitude);

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
    // r = a * (1 - e * e) / (1 + e);                  // équivalent
    // v = sqrt(STD_SUN / a / (1 - e * e)) * (1 + e);  // équivalent
    r = a * (1 - e);                                // est équivalent
    v = sqrt(STD_SUN * a * (1 - e * e)) / r;        // idem

    // State vectors
    body->pos.x =  r * cosw + sun.pos.x;
    body->pos.y =  r * sinw + sun.pos.y;
    body->vel.x = -v * sinw + sun.vel.x;
    body->vel.y =  v * cosw + sun.vel.y;
};