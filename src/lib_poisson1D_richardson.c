/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>  // Pour memcpy
#include <cblas.h>
#include <lapacke.h>
#include "lib_poisson1D.h"

void eig_poisson1D(double* eigval, int *la){
    for (int i = 1; i <= *la; i++) {
        eigval[i-1] = 4.0 * pow(sin((i*M_PI)/(2.0*(*la+1.0))), 2.0);
    }
}

double eigmax_poisson1D(int *la){
    double n = (double) *la;
    return 4.0 * pow(sin(M_PI/(2.0*(n+1.0))), 2.0);
}

double eigmin_poisson1D(int *la){
    double n = (double) *la;
    return 4.0 * pow(sin((n*M_PI)/(2.0*(n+1.0))), 2.0);
}

double richardson_alpha_opt(int *la){
    // α_opt = 2/(λmin + λmax)
    double lambda_min = eigmin_poisson1D(la);
    double lambda_max = eigmax_poisson1D(la);
    return 2.0/(lambda_min + lambda_max);
}

void richardson_alpha(double *AB, double *RHS, double *X, double *alpha, int *lab, int *la, int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite) {
    int incx = 1;
    double *AX = (double *)malloc(sizeof(double) * (*la));
    double *R = (double *)malloc(sizeof(double) * (*la));
    double norm_R0, norm_R;

    // Calcul de AX = AB*X (utilisation de BLAS)
    cblas_dgbmv(CblasColMajor, CblasNoTrans, 
                *la, *la, *kl, *ku, 
                1.0, AB, *lab, 
                X, incx, 
                0.0, AX, incx);

    // R = RHS - AX
    cblas_dcopy(*la, RHS, incx, R, incx);
    cblas_daxpy(*la, -1.0, AX, incx, R, incx);

    // Calcul de la norme du résidu initial
    norm_R0 = cblas_dnrm2(*la, R, incx);

    *nbite = 0;
    do {
        // AX = AB*X
        cblas_dgbmv(CblasColMajor, CblasNoTrans, 
                    *la, *la, *kl, *ku, 
                    1.0, AB, *lab, 
                    X, incx, 
                    0.0, AX, incx);

        // R = RHS - AX
        cblas_dcopy(*la, RHS, incx, R, incx);
        cblas_daxpy(*la, -1.0, AX, incx, R, incx);

        // X = X + alpha*R
        cblas_daxpy(*la, *alpha, R, incx, X, incx);

        // Calcul de la norme du résidu
        norm_R = cblas_dnrm2(*la, R, incx);
        resvec[*nbite] = norm_R/norm_R0;
        (*nbite)++;
    } while ((*nbite < *maxit) && (norm_R/norm_R0 > *tol));

    free(AX);
    free(R);
}

void jacobi_poisson1D(double *AB, double *RHS, double *X, int *lab, int *la, int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite) {
    int i, incx = 1;
    double *X_old = (double *)malloc(sizeof(double) * (*la));
    double *AX = (double *)malloc(sizeof(double) * (*la));
    double *R = (double *)malloc(sizeof(double) * (*la));
    double norm_R0, norm_R;

    // Calcul du résidu initial
    cblas_dgbmv(CblasColMajor, CblasNoTrans, 
                *la, *la, *kl, *ku, 
                1.0, AB, *lab, 
                X, incx, 
                0.0, AX, incx);
    
    cblas_dcopy(*la, RHS, incx, R, incx);
    cblas_daxpy(*la, -1.0, AX, incx, R, incx);
    norm_R0 = cblas_dnrm2(*la, R, incx);

    *nbite = 0;
    do {
        // Sauvegarde de l'ancienne solution
        cblas_dcopy(*la, X, incx, X_old, incx);

        // Itération de Jacobi avec BLAS
        for(i = 0; i < *la; i++) {
            double sum = RHS[i];
            if(i > 0) sum += X_old[i-1];
            if(i < *la-1) sum += X_old[i+1];
            X[i] = sum / 2.0;
        }

        // Calcul du nouveau résidu
        cblas_dgbmv(CblasColMajor, CblasNoTrans, 
                    *la, *la, *kl, *ku, 
                    1.0, AB, *lab, 
                    X, incx, 
                    0.0, AX, incx);
        
        cblas_dcopy(*la, RHS, incx, R, incx);
        cblas_daxpy(*la, -1.0, AX, incx, R, incx);
        norm_R = cblas_dnrm2(*la, R, incx);
        
        resvec[*nbite] = norm_R / norm_R0;
        (*nbite)++;
    } while ((*nbite < *maxit) && (norm_R / norm_R0 > *tol));

    free(X_old);
    free(AX);
    free(R);
}

void gauss_seidel_poisson1D(double *AB, double *RHS, double *X, int *lab, int *la, int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite) {
    int i, incx = 1;
    double *AX = (double *)malloc(sizeof(double) * (*la));
    double *R = (double *)malloc(sizeof(double) * (*la));
    double norm_R0, norm_R;

    // Calcul du résidu initial
    cblas_dgbmv(CblasColMajor, CblasNoTrans, 
                *la, *la, *kl, *ku, 
                1.0, AB, *lab, 
                X, incx, 
                0.0, AX, incx);
    
    cblas_dcopy(*la, RHS, incx, R, incx);
    cblas_daxpy(*la, -1.0, AX, incx, R, incx);
    norm_R0 = cblas_dnrm2(*la, R, incx);

    *nbite = 0;
    do {
        // Itération de Gauss-Seidel
        for(i = 0; i < *la; i++) {
            double sum = RHS[i];
            if(i > 0) sum += X[i-1];        
            if(i < *la-1) sum += X[i+1];    
            X[i] = sum / 2.0;
        }

        // Calcul du nouveau résidu
        cblas_dgbmv(CblasColMajor, CblasNoTrans, 
                    *la, *la, *kl, *ku, 
                    1.0, AB, *lab, 
                    X, incx, 
                    0.0, AX, incx);
        
        cblas_dcopy(*la, RHS, incx, R, incx);
        cblas_daxpy(*la, -1.0, AX, incx, R, incx);
        norm_R = cblas_dnrm2(*la, R, incx);
        
        resvec[*nbite] = norm_R / norm_R0;
        (*nbite)++;
    } while ((*nbite < *maxit) && (norm_R / norm_R0 > *tol));

    free(AX);
    free(R);
}

