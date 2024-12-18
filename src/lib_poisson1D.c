/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#define _USE_MATH_DEFINES
#include "lib_poisson1D.h"
#include <string.h> 
#include <time.h>    
#include <math.h>    

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
    for (int i = 0; i < (*lab) * (*la); i++) {
        AB[i] = 0.0;
    }
      // Sur-diagonale 
    for (int j = 0; j < *la - 1; j++) {
        AB[0 + j*(*lab)] = -1.0;
    }
    
    // Diagonale principale 
    for (int j = 0; j < *la; j++) {
        AB[1 + j*(*lab)] = 2.0;
    }   
    // Sous-diagonale 
    for (int j = 1; j < *la; j++) {
        AB[2 + j*(*lab)] = -1.0;
    }
}
void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
    int i;
    for (i = 0; i < (*lab) * (*la); i++) {
        AB[i] = 0.0;
    }
    
    for (i = 0; i < *la; i++) {
        AB[1 + i * (*lab)] = 1.0;
    }
}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1) {
    int i;
    double h = 1.0 / (*la + 1); 
    

    for (i = 0; i < *la; i++) {
        RHS[i] = h * h;  // f(x) * h^2
    }
    
    RHS[0] = RHS[0] + (*BC0);        
    RHS[*la - 1] = RHS[*la - 1] + (*BC1);  
}

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1) {
    int i;
    double h = 1.0 / (*la + 1);
    
    for (i = 0; i < *la; i++) {
        double x = X[i];
        EX_SOL[i] = -0.5 * x * (x - 1) + (*BC0) * (1 - x) + (*BC1) * x;
    }
}

void set_grid_points_1D(double* x, int* la) {
    int i;
    double h = 1.0 / (*la + 1);
    for (i = 0; i < *la; i++) {
        x[i] = (i + 1) * h;
    }
}

double relative_forward_error(double* x, double* y, int* la) {
    int i;
    double max_err = 0.0;
    double max_norm = 0.0;
    
    for (i = 0; i < *la; i++) {
        double err = fabs(x[i] - y[i]);
        double norm = fabs(y[i]);
        if (err > max_err) max_err = err;
        if (norm > max_norm) max_norm = norm;
    }
    
    if (max_norm < 1e-15) return max_err;  
    return max_err / max_norm;
}

int indexABCol(int i, int j, int *lab){
    return i + j * (*lab);
}

int dgbtrftridiag(int *la, int *n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
    if (*kl != 1 || *ku != 1) {
        *info = -1;
        return *info;
    }

    *info = 0;
    
    for(int j = 0; j < *n-1; j++) {
 
        if(fabs(AB[2 + j*(*lab)]) < DBL_EPSILON) {
            *info = j+1;
            return *info;
        }
        

        double multiplier = AB[3 + j*(*lab)] / AB[2 + j*(*lab)];
        
        // Mise à jour de la sous-diagonale
        AB[3 + j*(*lab)] = multiplier;
        
        // Mise à jour de la diagonale suivante
        AB[2 + (j+1)*(*lab)] -= multiplier * AB[1 + j*(*lab)];
    }
    

    if(fabs(AB[2 + (*n-1)*(*lab)]) < DBL_EPSILON) {
        *info = *n;
        return *info;
    }
    
    for(int i = 0; i < *n; i++) {
        ipiv[i] = i+1;
    }
    
    return *info;
}
void test_dgbmv_poisson1D(double* AB, int *lab, int *la, double *x, double *y) {
    enum CBLAS_ORDER order = CblasColMajor;
    enum CBLAS_TRANSPOSE trans = CblasNoTrans;
    int kl = 1;                // Une sous-diagonale
    int ku = 1;                // Une sur-diagonale
    double alpha = 1.0;
    double beta = 0.0;
    int incx = 1;
    int incy = 1;
    
    // Décalage pour pointer sur la première diagonale utile
    double* AB_offset = AB + 1;  // Pointer sur la sur-diagonale
    
    cblas_dgbmv(order, trans, 
                *la,        // M = nombre de lignes
                *la,        // N = nombre de colonnes
                kl,         // nombre de sous-diagonales
                ku,         // nombre de sur-diagonales
                alpha,
                AB_offset,  // matrice bande avec offset
                *lab,       // leading dimension de AB
                x,          // vecteur x
                incx,
                beta,
                y,          // vecteur résultat y
                incy);
}

void test_tridiag_LU(double* AB, int *lab, int *la) {
    printf("\n=== Test détaillé de la factorisation LU tridiagonale ===\n");
    
    int kl = 1;
    int ku = 1;
    int info;
    int NRHS = 1;
    char trans = 'N';
    
 
    double* AB_orig = (double*)malloc((*lab) * (*la) * sizeof(double));
    double* AB_lapack = (double*)malloc((*lab) * (*la) * sizeof(double));
    memcpy(AB_orig, AB, (*lab) * (*la) * sizeof(double));
    memcpy(AB_lapack, AB, (*lab) * (*la) * sizeof(double));
    
  
    double* x_exact = (double*)malloc((*la) * sizeof(double));
    double* b = (double*)malloc((*la) * sizeof(double));
    double* x_computed = (double*)malloc((*la) * sizeof(double));
    double* x_lapack = (double*)malloc((*la) * sizeof(double));
    int* ipiv = (int*)malloc((*la) * sizeof(int));
    int* ipiv_lapack = (int*)malloc((*la) * sizeof(int));
    
  
    for(int i = 0; i < *la; i++) {
        x_exact[i] = sin(i * 3.14159265358979323846 / (*la)) + 1.0;  // valeurs entre 0 et 2
    }

    test_dgbmv_poisson1D(AB_orig, lab, la, x_exact, b);
    
    printf("\n1. Test de la factorisation LU personnalisée :\n");
    printf("--------------------------------------\n");
    printf("Dimension de la matrice : %d x %d\n", *la, *la);
    printf("Largeur de bande : %d\n", *lab);
    

    clock_t start = clock();
    dgbtrftridiag(la, la, &kl, &ku, AB, lab, ipiv, &info);
    clock_t end = clock();
    double time_lu = ((double)(end - start)) / CLOCKS_PER_SEC;
    
    printf("Statut de la factorisation : %s (info = %d)\n", 
           (info == 0) ? "Succès" : "Échec", info);
    
    if(info != 0) {
        printf("ERREUR : Factorisation LU personnalisée échouée (info = %d)\n", info);
        goto cleanup;
    }
    
 
    memcpy(x_computed, b, (*la) * sizeof(double));
    start = clock();
    dgbtrs_(&trans, la, &kl, &ku, &NRHS, AB, lab, ipiv, x_computed, la, &info, 1);
    end = clock();
    double time_solve = ((double)(end - start)) / CLOCKS_PER_SEC;
    
  
    printf("\n2. Comparaison avec LAPACK :\n");
    printf("--------------------------------------\n");
    printf("Exécution de la factorisation LAPACK...\n");
    
    memcpy(x_lapack, b, (*la) * sizeof(double));
    start = clock();
    dgbtrf_(la, la, &kl, &ku, AB_lapack, lab, ipiv_lapack, &info);
    end = clock();
    double time_lu_lapack = ((double)(end - start)) / CLOCKS_PER_SEC;
    
    printf("Statut de la factorisation LAPACK : %s (info = %d)\n",
           (info == 0) ? "Succès" : "Échec", info);
    
    if(info == 0) {
        start = clock();
        dgbtrs_(&trans, la, &kl, &ku, &NRHS, AB_lapack, lab, ipiv_lapack, x_lapack, la, &info, 1);
        end = clock();
        double time_solve_lapack = ((double)(end - start)) / CLOCKS_PER_SEC;

        double max_error_custom = 0.0;
        double max_error_lapack = 0.0;
        double norm_x = 0.0;
        
        for(int i = 0; i < *la; i++) {
            max_error_custom = fmax(max_error_custom, fabs(x_computed[i] - x_exact[i]));
            max_error_lapack = fmax(max_error_lapack, fabs(x_lapack[i] - x_exact[i]));
            norm_x = fmax(norm_x, fabs(x_exact[i]));
        }
        
        printf("\n3. Résultats des tests :\n");
        printf("--------------------------------------\n");
        printf("Méthode personnalisée :\n");
        printf("- Temps factorisation : %e sec\n", time_lu);
        printf("- Temps résolution    : %e sec\n", time_solve);
        printf("- Erreur relative     : %e\n", max_error_custom/norm_x);
        
        printf("\nMéthode LAPACK :\n");
        printf("- Temps factorisation : %e sec\n", time_lu_lapack);
        printf("- Temps résolution    : %e sec\n", time_solve_lapack);
        printf("- Erreur relative     : %e\n", max_error_lapack/norm_x);
        
        // 7. Vérifier la structure
        printf("\n4. Vérification de la structure :\n");
        printf("--------------------------------------\n");
        int structure_ok = 1;
        for(int j = 0; j < *la; j++) {
            for(int i = 0; i < *lab; i++) {
                if(i > 3 || i < 1) {  // En dehors de la bande
                    if(fabs(AB[i + j * (*lab)]) > 1e-12) {
                        structure_ok = 0;
                        printf("Élément non nul hors bande : AB[%d,%d] = %e\n", i, j, AB[i + j * (*lab)]);
                    }
                }
            }
        }
        
        if(structure_ok) {
            printf("Structure tridiagonale préservée ✓\n");
        }
    } else {
        printf("ERREUR : Factorisation LAPACK échouée (info = %d)\n", info);
    }
    
cleanup:
    free(AB_orig);
    free(AB_lapack);
    free(x_exact);
    free(b);
    free(x_computed);
    free(x_lapack);
    free(ipiv);
    free(ipiv_lapack);
    
    printf("\n=== Fin du test ===\n\n");
}

