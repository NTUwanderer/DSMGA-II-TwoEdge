/***************************************************************************
 *   Copyright (C) 2015 Tian-Li Yu and Shih-Huan Hsu                       *
 *   tianliyu@ntu.edu.tw                                                   *
 ***************************************************************************/


#include <math.h>
#include <iostream>
#include <cstdlib>
#include <time.h>
#include "statistics.h"
#include "dsmga2.h"
#include "global.h"
#include "chromosome.h"

using namespace std;


int
main (int argc, char *argv[]) {
    if (argc < 9) {
        printf ("DSMGA2 ell nInitial function maxGen maxFe repeat display rand_seed s_num=1 nk_step=1\n");
        printf ("function: \n");
        printf ("     ONEMAX:  0\n");
        printf ("     MK    :  1\n");
        printf ("     FTRAP :  2\n");
        printf ("     CYC   :  3\n");
        printf ("     NK    :  4\n");
        printf ("     SPIN  :  5\n");
        printf ("     SAT   :  6\n");
        printf ("     L_FTRAP  :  7\n");
        printf ("     L_2FTRAP :  8\n");
        printf ("     MKP   :  9\n");

        return -1;
    }

    int ell = atoi (argv[1]); // problem size
    int nInitial = atoi (argv[2]); // initial population size
    int fffff = atoi (argv[3]); // function
    int maxGen = atoi (argv[4]); // max generation
    int maxFe = atoi (argv[5]); // max fe
    int repeat = atoi (argv[6]); // how many time to repeat
    int display = atoi (argv[7]); // display each generation or not
    int rand_seed = atoi (argv[8]);  // rand seed
	int s_num = argc > 9 ? atoi (argv[9]) : 1;
	int nk_step = argc > 10 ? atoi (argv[10]) : 1;


    if (fffff == 4) {

        char filename[200];
        //sprintf(filename, "./NK_Instance/pnk%d_%d_%d_%d", ell, 4, 1, 1);
        sprintf(filename, "./NK_Instance/pnk%d_%d_%d_%d", ell, 4, nk_step, s_num);

        if (SHOW_BISECTION) printf("Loading: %s\n", filename);
        FILE *fp = fopen(filename, "r");
        loadNKWAProblem(fp, &nkwa);
        fclose(fp);
    }

    if (fffff == 5) {
        char filename[200];
        sprintf(filename, "./SPIN/%d/%d_%d",ell, ell, s_num);
        if (SHOW_BISECTION) printf("Loading: %s\n", filename);
        loadSPIN(filename, &mySpinGlassParams);
    }

    if (fffff == 6) {
        char filename[200];
        sprintf(filename, "./SAT/uf%d/uf%d-0%d.cnf", ell, ell, s_num);
        if (SHOW_BISECTION) printf("Loading: %s\n", filename);
        loadSAT(filename, &mySAT);
    }

    if (fffff == 9) {
        char filename[200];
        sprintf(filename, "./All-MKP-Instances/sac94/weish/weish%02d.dat", s_num);
        if (SHOW_BISECTION) printf("Loading: %s\n", filename);
        loadMKP(filename, &myMKP);
        ell = myMKP.var;
    }


    if (rand_seed != -1)  // time
        myRand.seed((unsigned long)rand_seed);

    int i;

    Statistics stGen, stFE, stLSFE;
    int usedGen;

    int failNum = 0;
    double* bestFs = new double[repeat];

    for (i = 0; i < repeat; i++) {

        DSMGA2 ga (ell, nInitial, maxGen, maxFe, fffff);

        if (display == 1)
            usedGen = ga.doIt (true);
        else
            usedGen = ga.doIt (false);

        bestFs[i] = ga.bestF();

        if (!ga.foundOptima()) {
            failNum++;
            printf ("-");
        } else {
            stFE.record (Chromosome::hitnfe);
            stLSFE.record (Chromosome::lsnfe);
            stGen.record (usedGen);
            printf ("+");
        }

        fflush (NULL);

    }
    cout<<endl; 
    // printf ("%f  %f  %f %d\n", stGen.getMean (), stFE.getMean(), stLSFE.getMean(), failNum);
    printf ("\n");
    printf ("Gen: %f\n", stGen.getMean ());
    printf ("FailNum: %d\n", failNum);
    printf ("LSNFE: %f\n", stLSFE.getMean());
    printf ("NFE: %f\n", stFE.getMean());
    printf ("N_std: %f\n", stFE.getStdev());
    printf ("Fs:");
    for (int i = 0; i < repeat; ++i)
        printf (" %f", bestFs[i]);
    printf ("\n");

    delete[] bestFs;

    if (fffff == 4) freeNKWAProblem(&nkwa);
    if (fffff == 9) freeMKPinstance(&myMKP);

    return EXIT_SUCCESS;
}
