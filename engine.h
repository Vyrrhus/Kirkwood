#ifndef ENGINE_H
    #define ENGINE_H
    
    #include "kirkwood.h"

    /***** PROTOTYPES */
    // Integrators
    Body   leapfrog(int nbMajorBodies, Body listMajorBodies[], Body listAsteroids[], double dT);
    Body    yoshida(int nbMajorBodies, Body listMajorBodies[], Body listAsteroids[], double dT);
    Body yoshida8th(int nbMajorBodies, Body listMajorBodies[], Body listAsteroids[], double dT);

    // Physic model
    void kepler_from_state(Body sun, Body *body);
    void state_from_kepler(Body sun, Body *body);

#endif