/******************************************/
/* tp2_poisson1D_direct.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "lib_poisson1D.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define TRF 0
#define TRI 1
#define SV 2

int main(int argc, char *argv[])
{
    int nbpoints = 10;
    int la = nbpoints - 2;
    int lab = 4;
    int kv = 1;
    int NRHS = 1;
    char trans = 'N';
    double T0 = -5.0;
    double T1 = 5.0;

    printf("--------- Poisson 1D ---------\n\n");
    printf("Configuration :\n");
    printf("- Taille matrice : %d x %d\n", la, la);
    printf("- Largeur bande : %d\n", lab);
    printf("- Conditions : T0 = %f, T1 = %f\n\n", T0, T1);

    // Allocation et initialisation
    double *AB = (double *)malloc(lab * la * sizeof(double));
    double *RHS = (double *)malloc(la * sizeof(double));
    int *ipiv = (int *)malloc(la * sizeof(int));
    
    // Construction du système
    set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
    set_dense_RHS_DBC_1D(RHS, &la, &T0, &T1);

    // === Exercice 4 : Test DGBMV ===
    printf("=== Exercice 4 : Test DGBMV ===\n");
    double *x = (double *)malloc(la * sizeof(double));
    double *y = (double *)malloc(la * sizeof(double));
    for(int i = 0; i < la; i++) x[i] = 1.0;
    
    clock_t start = clock();
    test_dgbmv_poisson1D(AB, &lab, &la, x, y);
    clock_t end = clock();
    
    printf("1. Validation :\n");
    printf("y = Ax avec x = [1,...,1] doit donner [1,...,1]\n");
    for(int i = 0; i < la; i++) {
        printf("y[%d] = %f\n", i, y[i]);
    }
    printf("\n2. Performance DGBMV :\n");
    printf("Temps : %e sec\n", ((double)(end-start))/CLOCKS_PER_SEC);
    printf("Complexité théorique : O(n) pour matrice bande\n\n");



    // === Exercices 5-6 : Méthodes directes ===
    printf("=== Exercice 5-6 : Méthodes directes ===\n");
    if(argc > 1) {
        int method = atoi(argv[1]);

        printf("=== Test méthode %d ===\n", method);
        printf("0 = DGBTRF+DGBTRS (Ex.5)\n\n");
        printf("1 = LU tridiagonal (Ex.6)\n");
        printf("2 = DGBSV (Ex.5)\n");


        int info = 0;
        start = clock();
        
        switch(method) {
            case TRI:  // 1
                printf("Méthode : LU tridiagonal\n");
                printf("Complexité théorique : O(n)\n");
                dgbtrftridiag(&la, &la, &kv, &kv, AB, &lab, ipiv, &info);
                break;
            case SV:   // 2
                printf("Méthode : DGBSV\n");
                printf("Complexité théorique : O(n) pour matrice bande\n");
                dgbsv_(&la, &kv, &kv, &NRHS, AB, &lab, ipiv, RHS, &la, &info);
                break;
            case TRF:  // 0
                printf("Méthode : DGBTRF + DGBTRS\n");
                printf("Complexité théorique : O(n) pour matrice bande\n");
                dgbtrf_(&la, &la, &kv, &kv, AB, &lab, ipiv, &info);
                if(info == 0) {
                    char trans = 'N';
                    int ldab = lab;
                    int ldb = la;
                    dgbtrs_(&trans, &la, &kv, &kv, &NRHS, AB, &ldab, ipiv, RHS, &ldb, &info, 1);
                }
                break;
        }
        end = clock();
        
        printf("\nRésultats :\n");
        printf("- Statut : %s (info = %d)\n", (info == 0) ? "Succès" : "Échec", info);
        printf("- Temps total : %e sec\n", ((double)(end-start))/CLOCKS_PER_SEC);
    }

    free(AB); free(RHS); free(ipiv); free(x); free(y);
    printf("\n--------- Fin -----------\n");
    return EXIT_SUCCESS;
}
