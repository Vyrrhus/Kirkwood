#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "kirkwood.h"
#include "engine.h"
#include "file.h"

/***** PROTOTYPES */
static Body run_simulation(double *time, int nbMajorBodies, Body listMajorBodies[], Body listAsteroids[]);

/***** MAIN */
int main() {
    double time = 0;

    // Generate asteroid input file
    if (RESET_ASTEROID_INPUT) {
        generate_asteroid();
    }

    // Initiate Body tabs
    double nbMajorBodies = get_nb_major_body();
    Body asteroids[(int) fmax(ASTEROID_NUMBER, 1)];     // Minimum to 1 to avoid errors
    Body* majorBodies;
    Body  sun;

    majorBodies = (Body*) malloc(nbMajorBodies * sizeof(Body));
    if (majorBodies == NULL) {
        printf("Dynamic allocation failed !");
        return 1;
    }
    sun = init_input(nbMajorBodies, majorBodies, asteroids);

    // Initiate output file
    if (SAVE_BOOL) {
        init_outputFile(nbMajorBodies);
        save(time, nbMajorBodies, majorBodies, asteroids, sun);
    }

    // Run simulation
    sun = run_simulation(&time, nbMajorBodies, majorBodies, asteroids);

    // Save initial and final Kepler elements for the asteroids
    if (SAVE_HIST) {
        save_hist(nbMajorBodies, sun, asteroids);
    }
    
    return 0;
}

/***** FUNCTIONS */
static Body run_simulation(double *time, int nbMajorBodies, Body listMajorBodies[], Body listAsteroids[]) {
    /*
     *  Integrate the system of N bodies
     *  Length time: [0 - LENGTH_SIMULATION]
     *  Step time  : STEP_TIME
     *  
     *  Integrator : Leapfrog / Yoshida
    */
   int stop = 0;
   double percent;
   Body sun;

   printf("BEGIN SIMULATION\n");
   while (!stop) {
       *time += STEP_TIME;

       if (*time >= LENGTH_SIMULATION) {
           *time = LENGTH_SIMULATION;
           stop  = 1;
       }

       // Iteration
    //    sun = leapfrog(nbMajorBodies, listMajorBodies, listAsteroids, (double) STEP_TIME);
       sun = yoshida(nbMajorBodies, listMajorBodies, listAsteroids, (double) STEP_TIME);

       // Save data
       if (fmod(*time, (double) fmax(STEP_TIME,STEP_WRITING)) < STEP_TIME || *time == LENGTH_SIMULATION) {
           save(*time, nbMajorBodies, listMajorBodies, listAsteroids, sun);
           percent = *time / ((double) LENGTH_SIMULATION) * 100.;
           printf("\rSimulation: %.5f%%", percent);
       }
   }

   printf("\nEND OF SIMULATION");

   return sun;
};