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
    Vector acceleration = {0,0,0};
    Vector relativePosition;
    double r2;
    double a;

    for (int i = 0 ; i < nbMajorBodies ; i++) {
        relativePosition.x = body.pos.x - listMajorBodies[i].pos.x;
        relativePosition.y = body.pos.y - listMajorBodies[i].pos.y;
        relativePosition.z = body.pos.z - listMajorBodies[i].pos.z;

        r2 = relativePosition.x * relativePosition.x + relativePosition.y * relativePosition.y + relativePosition.z * relativePosition.z;
        a  = listMajorBodies[i].std / (r2 * sqrt(r2));

        if (r2 != 0) {
            acceleration.x -= relativePosition.x * a;
            acceleration.y -= relativePosition.y * a;
            acceleration.z -= relativePosition.z * a;
        }
    }

    return acceleration;
};

// Integrators
Body leapfrog(int nbMajorBodies, Body listMajorBodies[], Body listAsteroids[], double dT) {
    /*
     *  Leapfrog integrator
     */
    Vector acc = {0, 0, 0};
    Body sun;

    // Position at t+dt
    for (int i = 0 ; i < ASTEROID_NUMBER ; i++) {
        listAsteroids[i].pos.x += listAsteroids[i].vel.x * dT + 0.5 * listAsteroids[i].acc.x * dT * dT;
        listAsteroids[i].pos.y += listAsteroids[i].vel.y * dT + 0.5 * listAsteroids[i].acc.y * dT * dT;
        listAsteroids[i].pos.z += listAsteroids[i].vel.z * dT + 0.5 * listAsteroids[i].acc.z * dT * dT;
    }
    for (int i = 0 ; i < nbMajorBodies ; i++) {
        listMajorBodies[i].pos.x += listMajorBodies[i].vel.x * dT + 0.5 * listMajorBodies[i].acc.x * dT * dT;
        listMajorBodies[i].pos.y += listMajorBodies[i].vel.y * dT + 0.5 * listMajorBodies[i].acc.y * dT * dT;
        listMajorBodies[i].pos.z += listMajorBodies[i].vel.z * dT + 0.5 * listMajorBodies[i].acc.z * dT * dT;
    }

    // Velocity and acceleration at t+dt
    for (int i = 0 ; i < ASTEROID_NUMBER ; i++) {
        acc = get_acceleration(listAsteroids[i], nbMajorBodies, listMajorBodies);
        listAsteroids[i].vel.x += 0.5 * (listAsteroids[i].acc.x + acc.x) * dT;
        listAsteroids[i].vel.y += 0.5 * (listAsteroids[i].acc.y + acc.y) * dT;
        listAsteroids[i].vel.z += 0.5 * (listAsteroids[i].acc.z + acc.z) * dT;
        listAsteroids[i].acc = acc;
    }
    for (int i = 0 ; i < nbMajorBodies ; i++) {
        acc = get_acceleration(listMajorBodies[i], nbMajorBodies, listMajorBodies);
        listMajorBodies[i].vel.x += 0.5 * (listMajorBodies[i].acc.x + acc.x) * dT;
        listMajorBodies[i].vel.y += 0.5 * (listMajorBodies[i].acc.y + acc.y) * dT;
        listMajorBodies[i].vel.z += 0.5 * (listMajorBodies[i].acc.z + acc.z) * dT;
        listMajorBodies[i].acc = acc;

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
    Vector acc = {0, 0, 0};

    // Constants for Yoshida steps
    double a = pow(2, (double) 1/3);
    double W0 = - a / (2 - a);
    double W1 =   1 / (2 - a);
    double c[4] = {W1/2, (W0+W1)/2, (W0+W1)/2, W1/2};
    double d[3] = {W1, W0, W1};

    // Yoshida steps
    for (int step = 0 ; step < 4 ; step++) {
        // Position
        for (int i = 0 ; i < ASTEROID_NUMBER ; i++) {
            if (listAsteroids[i].semiMajorAxis >= SMA_MAX_VALUE || listAsteroids[i].eccentricity >= ECC_MAX_VALUE) {
                continue;
            }
            listAsteroids[i].pos.x += c[step] * listAsteroids[i].vel.x * dT;
            listAsteroids[i].pos.y += c[step] * listAsteroids[i].vel.y * dT;
            listAsteroids[i].pos.z += c[step] * listAsteroids[i].vel.z * dT;
        }
        for (int i = 0 ; i < nbMajorBodies ; i++) {
            listMajorBodies[i].pos.x += c[step] * listMajorBodies[i].vel.x * dT;
            listMajorBodies[i].pos.y += c[step] * listMajorBodies[i].vel.y * dT;
            listMajorBodies[i].pos.z += c[step] * listMajorBodies[i].vel.z * dT;
        }

        // No need to compute velocity for the 4th step
        if (step == 3) {
            break;
        }

        // Velocity and acceleration
        for (int i = 0 ; i < ASTEROID_NUMBER ; i++) {
            if (listAsteroids[i].semiMajorAxis > SMA_MAX_VALUE || listAsteroids[i].eccentricity >= ECC_MAX_VALUE) {
                continue;
            }
            acc = get_acceleration(listAsteroids[i], nbMajorBodies, listMajorBodies);
            listAsteroids[i].vel.x += d[step] * acc.x * dT;
            listAsteroids[i].vel.y += d[step] * acc.y * dT;
            listAsteroids[i].vel.z += d[step] * acc.z * dT;
            listAsteroids[i].acc = acc;
        }
        for (int i = 0 ; i < nbMajorBodies ; i++) {
            acc = get_acceleration(listMajorBodies[i], nbMajorBodies, listMajorBodies);
            listMajorBodies[i].vel.x += d[step] * acc.x * dT;
            listMajorBodies[i].vel.y += d[step] * acc.y * dT;
            listMajorBodies[i].vel.z += d[step] * acc.z * dT;
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
    double sunStdParameter = (double) STD_SUN;

    // Relative state vectors to the Sun
    Vector relPos = {body->pos.x - sun.pos.x, 
                     body->pos.y - sun.pos.y,
                     body->pos.z - sun.pos.z};
    Vector relVel = {body->vel.x - sun.vel.x, 
                     body->vel.y - sun.vel.y,
                     body->vel.z - sun.vel.z};

    // Squared distance and velocity
    double r2 = relPos.x * relPos.x + relPos.y * relPos.y + relPos.z * relPos.z;
    double v2 = relVel.x * relVel.x + relVel.y * relVel.y + relVel.z * relVel.z;

    // Sun
    if (r2 == 0) {
        return;
    }
    double rNorm = sqrt(r2);

    // Orbital momentum vector : h = r ^ v
    Vector h = {relPos.y * relVel.z - relPos.z * relVel.y,
                relPos.z * relVel.x - relPos.x * relVel.z,
                relPos.x * relVel.y - relPos.y * relVel.x};
    double hNorm = sqrt(h.x * h.x + h.y * h.y + h.z * h.z);

    // Node vector : n = K ^ h
    Vector n = { -h.y,
                  h.x,
                  0};
    double nNorm = sqrt(n.x * n.x + n.y * n.y + n.z * n.z);

    // Eccentricity vector : e = v ^ h / std - r / |r|
    Vector e = {(relVel.y * h.z - relVel.z * h.y) / sunStdParameter - relPos.x / rNorm,
                (relVel.z * h.x - relVel.x * h.z) / sunStdParameter - relPos.y / rNorm,
                (relVel.x * h.y - relVel.y * h.x) / sunStdParameter - relPos.z / rNorm};
    double eNorm = sqrt(e.x * e.x + e.y * e.y + e.z * e.z);
    
    // Mechanical specific energy & semi-major axis (given eNorm != 1)
    double energy = v2 / 2 - sunStdParameter / rNorm;
    double semiMA = - sunStdParameter / (2 * energy);

    // Angles
    double i = acos(h.z / hNorm);
    double W = acos(n.x / nNorm);
    if (n.y < 0) {
        W = 2 * M_PI - W;
    }
    double w = acos((e.x * relPos.x + e.y * relPos.y + e.z * relPos.z) / (eNorm * rNorm));

    // Special angles
    double wTrue = acos(e.x / eNorm);
    if (e.y < 0) {
        wTrue = 2 * M_PI - wTrue;
    }

    // Elliptical equatorial (ie i = 0)
    if (i <= 1e-12 && eNorm > 1e-12) {
       W = 0;
       w = wTrue;
    } 
    
    // Circular inclined (ie e = 0)
    if (i > 1e-12 && eNorm <= 1e-12) {
        w = 0;
    }

    // Circular equatorial (e = 0, i = 0)
    if ( i <= 1e-12 && eNorm <= 1e-12) {
        w = 0;
        W = 0;
    }

    body->semiMajorAxis = semiMA;
    body->eccentricity  = eNorm;
    body->inclination   = i;
    body->longitudeAN   = W;
    body->argPeriapsis  = w;
};
void state_from_kepler(Body sun, Body *body) {
    /*  Compute state vectors (pos, vel) from Kepler elements (a, e, i, W, w).
        We assume the true anomaly is 0 (ie the body is defined at its periapsis)
     */
    double a, e, i, W, w;
    double cosw, sinw, cosi, sini, cosW, sinW;
    double std = (double) STD_SUN;

    a = body->semiMajorAxis;
    e = body->eccentricity;
    i = body->inclination;
    W = body->longitudeAN;
    w = body->argPeriapsis;

    cosi = cos(i);
    sini = sin(i);
    cosW = cos(W);
    sinW = sin(W);
    cosw = cos(w);
    sinw = sin(w);

    // Distance & velocity relative to the Sun
    double r = a * (1 - e);
    double v = sqrt(std / (a * (1 - e*e))) * (1 + e);

    // State vectors
    Vector pos = {r * (cosW*cosw - sinW*sinw*cosi),
                  r * (sinW*cosw + cosW*sinw*cosi),
                  r * (sinw * sini)};
    
    Vector vel = {v * (-cosW*sinw - sinW*cosw*cosi),
                  v * (-sinW*sinw + cosW*cosw*cosi),
                  v * (cosw * sini)};

    // Relative state vectors
    body->pos = pos;
    body->vel = vel;

    // True state vectors
    body->pos.x += sun.pos.x;
    body->pos.y += sun.pos.y;
    body->pos.z += sun.pos.z;

    body->vel.x += sun.vel.x;
    body->vel.y += sun.vel.y;
    body->vel.z += sun.vel.z;
};