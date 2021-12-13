#ifndef FILE_H
    #define FILE_H

    #include "kirkwood.h"

    /***** PROTOTYPES */
    // Output
    void init_outputFile(int nbMajorBodies);
    void save(double time, int nbMajorBodies, Body listMajorBodies[], Body listAsteroid[], Body sun);
    void save_hist(int nbMajorBodies, Body sun, Body listAsteroid[]);

    // Input
    void generate_asteroid();
    int  get_nb_major_body();
    Body init_input(int nbMajorBodies, Body listMajorBodies[], Body listAsteroid[]);

#endif