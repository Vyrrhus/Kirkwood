#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "kirkwood.h"
#include "file.h"
#include "engine.h"

/***** FILE NAME */
static char* INPUT_BODIES     = "input/saturn.dat";
static char* INPUT_ASTEROID   = "input/asteroids.dat";
static char* OUTPUT_FILENAME  = "output3d/%dMB_%drev_%dast.dat";
static char* HIST_FILENAME    = "output3d/%dMB_%drev_%dast-hist.dat";

/***** FUNCTIONS */
// Static functions
static double randfrom(double min, double max) {
    /*  Return random number between [min, max]
     */
    double range = (max - min);
    double div   = RAND_MAX / range;
    return min + (rand() / div);
};
static void init_asteroid(Body sun, Body listAsteroid[]) {
    FILE *datafile = NULL;
    int id;
    Vector null = {0, 0};

    datafile = fopen(INPUT_ASTEROID, "r");
    if (datafile == NULL) {
        printf("Error while reading INPUT_ASTEROID file");
        exit(EXIT_FAILURE);
    }

    // printf("ASTEROID INIT\n");
    for (int i = 0 ; i < ASTEROID_NUMBER ; i++) {
        fscanf(datafile, "%d,%lf,%lf,%lf\n",
               &id,
               &listAsteroid[i].semiMajorAxis,
               &listAsteroid[i].eccentricity,
               &listAsteroid[i].trueLongitude);
        state_from_kepler(sun, &listAsteroid[i]);
        listAsteroid[i].acc = null;

    }
};
static Body init_major_bodies(int nbMajorBodies, Body listMajorBodies[]) {
    FILE *datafile = NULL;
    char name[50];
    Body sun;
    Vector null = {0, 0};

    datafile = fopen(INPUT_BODIES, "r");
    if (datafile == NULL) {
        printf("Error while reading INPUT_BODIES file");
        exit(EXIT_FAILURE);
    }

    // Sun is always the first line with format string:
    // {name,std,x,y,z,vx,vy,vz}
    fscanf(datafile, "%[^,],%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",
           &name,&listMajorBodies[0].std,
           &listMajorBodies[0].pos.x,&listMajorBodies[0].pos.y,&listMajorBodies[0].pos.z,
           &listMajorBodies[0].vel.x,&listMajorBodies[0].vel.y,&listMajorBodies[0].vel.z);
    
    listMajorBodies[0].semiMajorAxis = 0;
    listMajorBodies[0].eccentricity  = 0;
    listMajorBodies[0].inclination   = 0;
    listMajorBodies[0].longitudeAN   = 0;
    listMajorBodies[0].argPeriapsis  = 0;
    listMajorBodies[0].std          *= GM_CONVERSION;
    listMajorBodies[0].acc           = null;
        }

    // Check if Sun is well declared
    sun = listMajorBodies[0];
    if (sun.std != STD_SUN) {
        exit(EXIT_FAILURE);
    }

    // Other major bodies follow with format string:
    // {name, std, a, e, i, W, w}
    for (int i = 1 ; i < nbMajorBodies ; i++) {
        fscanf(datafile, "%[^,],%lf,%lf,%lf,%lf,%lf,%lf\n",
               &name,
               &listMajorBodies[i].std,
               &listMajorBodies[i].semiMajorAxis,
               &listMajorBodies[i].eccentricity,
               &listMajorBodies[i].inclination,
               &listMajorBodies[i].longitudeAN,
               &listMajorBodies[i].argPeriapsis);

        listMajorBodies[i].inclination  *= RAD_CONVERSION;
        listMajorBodies[i].longitudeAN  *= RAD_CONVERSION;
        listMajorBodies[i].argPeriapsis *= RAD_CONVERSION;
        listMajorBodies[i].std          *= GM_CONVERSION;
        listMajorBodies[i].acc           = null;

        state_from_kepler(sun, &listMajorBodies[i]);

        }
    }
    return sun;
};

// Output
void init_outputFile(int nbMajorBodies) {
    /*
     *  Set the OUTPUT files name and reset them if they already exist
     */
    FILE *datafile = NULL;

    int length = snprintf(NULL, 0, OUTPUT_FILENAME, nbMajorBodies, JUPITER_REV, ASTEROID_NUMBER);
    char* str  = malloc(length+1);
    snprintf(str, length+1, OUTPUT_FILENAME, nbMajorBodies, JUPITER_REV, ASTEROID_NUMBER);
    OUTPUT_FILENAME = str;

    datafile = fopen(OUTPUT_FILENAME, "w");
    if (datafile == NULL) {
        printf("Error while opening OUTPUT file");
        exit(EXIT_FAILURE);
    }

    // Set header
    fprintf(datafile, "%d,id,x,y,vx,vy,semi-major axis,eccentricity,true longitude", ASTEROID_NUMBER + nbMajorBodies);
    // fprintf(datafile, "%d,id,semi-major axis,eccentricity,true longitude", ASTEROID_NUMBER + nbMajorBodies);

    fclose(datafile);
};
void save(double time, int nbMajorBodies, Body listMajorBodies[], Body listAsteroid[], Body sun) {
    /*
     *  Save data in OUTPUT file
     */
    if (!SAVE_BOOL) {
        return;
    }
    FILE *datafile = NULL;

    datafile = fopen(OUTPUT_FILENAME, "a");
    if (datafile == NULL) {
        printf("Error while opening OUTPUT file.");
        exit(EXIT_FAILURE);
    }

    // New line with time
    fprintf(datafile, "\n%.5e", time);

    // Add data for each major Body
    for (int i = 0 ; i < nbMajorBodies ; i++) {
        kepler_from_state(sun, &listMajorBodies[i]);
        fprintf(datafile, ",%d,%.5e,%.5e,%.5e,%.5e,%.5e,%.5e,%.5e", -(i+1),
                          listMajorBodies[i].pos.x, listMajorBodies[i].pos.y,
                          listMajorBodies[i].vel.x, listMajorBodies[i].vel.y,
                          listMajorBodies[i].semiMajorAxis, listMajorBodies[i].eccentricity, listMajorBodies[i].trueLongitude);
        // fprintf(datafile, ",%d,%.5e,%.5e,%.5e", -(i+1), listMajorBodies[i].semiMajorAxis,
        //                                                 listMajorBodies[i].eccentricity,
        //                                                 listMajorBodies[i].trueLongitude);
    }

    // Add data for each asteroid
    for (int i = 0 ; i < ASTEROID_NUMBER ; i++) {
        kepler_from_state(sun, &listAsteroid[i]);

        fprintf(datafile, ",%d,%.5e,%.5e,%.5e,%.5e,%.5e,%.5e,%.5e", (i+1),
                          listAsteroid[i].pos.x, listAsteroid[i].pos.y,
                          listAsteroid[i].vel.x, listAsteroid[i].vel.y,
                          listAsteroid[i].semiMajorAxis, listAsteroid[i].eccentricity, listAsteroid[i].trueLongitude);
        // fprintf(datafile, ",%d,%.5e,%.5e,%.5e", (i+1), listAsteroid[i].semiMajorAxis,
        //                                                listAsteroid[i].eccentricity,
        //                                                listAsteroid[i].trueLongitude);
    }

    fclose(datafile);

};
void save_hist(int nbMajorBodies, Body sun, Body listAsteroid[]) {
    /*
     *  Save asteroids elements at the end of the simulation
     */
    FILE *outputFile = NULL;
    FILE *inputFile  = NULL;
    int id;
    double a,e,w;

    int length = snprintf(NULL, 0, HIST_FILENAME, nbMajorBodies, JUPITER_REV, ASTEROID_NUMBER);
    char* str  = malloc(length+1);
    snprintf(str, length+1, HIST_FILENAME, nbMajorBodies, JUPITER_REV, ASTEROID_NUMBER);
    HIST_FILENAME = str;

    outputFile = fopen(HIST_FILENAME, "w");
    inputFile  = fopen(INPUT_ASTEROID,"r");
    if (outputFile == NULL || inputFile == NULL) {
        printf("Error while opening file");
        exit(EXIT_FAILURE);
    }

    int firstLine = 1;
    for (int i = 0 ; i < ASTEROID_NUMBER ; i++) {
        fscanf(inputFile, "%d,%lf,%lf,%lf\n",&id,&a,&e,&w);

        if (listAsteroid[i].isTooFar) {
            continue;
        }
        if (firstLine) {
            firstLine = 0;
        } else {
            fprintf(outputFile,"\n");
        }
        kepler_from_state(sun, &listAsteroid[i]);
        fprintf(outputFile,"%.5e,%.5e,%.5e,%.5e,%.5e,%.5e",
                            a,e,w,  listAsteroid[i].semiMajorAxis,
                                    listAsteroid[i].eccentricity,
                                    listAsteroid[i].trueLongitude);
    }
    fclose(outputFile);
    fclose(inputFile);
};

// Input
void generate_asteroid() {
    /*
     *  Generate datafile with {ASTEROID_NUMBER} lines :
     *  {id, semi-major axis, eccentricity, true longitude of periapsis}
     */
    FILE *datafile = NULL;
    double semiMajorAxis, eccentricity, trueLongitude;
    double r1, r2, peri, apo;

    datafile = fopen(INPUT_ASTEROID, "w");
    if (datafile == NULL) {
        printf("Error while opening INPUT_ASTEROID file");
        exit(EXIT_FAILURE);
    }

    // Initiate seed for rand()
    srand(time (NULL));

    for (int i = 1 ; i <= ASTEROID_NUMBER ; i++) {
        /*  Generate data
         *      Uniform distribution for semi-major axis within the given range
         *      Eccentricity set to 0
         *      Random true longitude for periapsis between 0 and 2*Pi
         */
        // semiMajorAxis = randfrom(MIN_AST_GENERATION, MAX_AST_GENERATION);
        // eccentricity  = 0;
        // semiMajorAxis = (double) (i-1) / (double) (ASTEROID_NUMBER - 1) * (MAX_AST_GENERATION - MIN_AST_GENERATION) + MIN_AST_GENERATION;
        r1 = randfrom(1.666, 4.951);
        r2 = randfrom(1.666, 4.951);
        if (r1 < r2) {
            peri = r1;
            apo  = r2;
        } else {
            peri = r2;
            apo  = r1;
        }
        eccentricity  = (1 - peri/apo) / (1 + peri/apo);
        semiMajorAxis = peri / (1 - eccentricity);
        trueLongitude = randfrom(0., 2*M_PI);
        
        // Save data
        fprintf(datafile, "%d, %.15e, %.15e, %.15e", i, semiMajorAxis, eccentricity, trueLongitude);
        if (i != ASTEROID_NUMBER) {
            fprintf(datafile, "\n");
        }
    }
    printf("%d asteroids generated.\n", ASTEROID_NUMBER);

    fclose(datafile);
};
int  get_nb_major_body() {
    /*
     *  Return number of major bodies as defined in INPUT_BODIES file
     *  Useful for dynamic allocation
     */
    FILE *datafile = NULL;
    int nbBody = 0;
    
    datafile = fopen(INPUT_BODIES, "r");
    if (datafile == NULL) {
        printf("Error while opening INPUT_BODIES file");
        exit(EXIT_FAILURE);
    }

    int ch = 0;
    while (!feof(datafile)) {
        ch = fgetc(datafile);
        if (ch == '\n' || feof(datafile)) {
            nbBody++;
        }
    }

    fclose(datafile);

    return nbBody;
};
Body init_input(int nbMajorBodies, Body listMajorBodies[], Body listAsteroid[]) {
    /*
     *  Initialize all the tabs from the input files
     */
    Body sun;
    sun = init_major_bodies(nbMajorBodies, listMajorBodies);
    init_asteroid(sun, listAsteroid);

    return sun;
}