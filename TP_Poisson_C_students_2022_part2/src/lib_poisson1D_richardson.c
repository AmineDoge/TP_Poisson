/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void eig_poisson1D(double* eigval, int *la){
    *eigval = 2.0/(eigmin_poisson1D(la) + eigmax_poisson1D(la));
}

double eigmax_poisson1D(int *la){
  return -2.0*cos((*la)*M_PI/(*la+1)) + 2.0 ;
}

double eigmin_poisson1D(int *la){
  return -2.0*cos(M_PI/(*la+1)) + 2.0 ;
}

double richardson_alpha_opt(int *la) {
    double eig_max = eigmax_poisson1D(la);
    double eig_min = eigmin_poisson1D(la);

    if (eig_max > 0 && eig_min > 0) {
        return 2.0 / (eig_max + eig_min);
    } else {
        fprintf(stderr, "Erreur : valeurs propres invalides\n");
        return 0.0;
    }
}


void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
    // Allocation de vecteurs temporaires
    double *resid = (double *)malloc(sizeof(double) * (*la));
    double norm_res, norm_rhs;

    // Initialisation
    *nbite = 0;

    // Calcul de la norme initiale de RHS (||b||)
    norm_rhs = cblas_dnrm2(*la, RHS, 1);

    // Initialisation du résidu : r = b - A*x
    cblas_dcopy(*la, RHS, 1, resid, 1); // resid = RHS
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, resid, 1);

    // Calcul de la norme initiale du résidu
    norm_res = cblas_dnrm2(*la, resid, 1) / norm_rhs;
    resvec[0] = norm_res;

    // Itérations de Richardson
    while (*nbite < *maxit && norm_res > *tol) {
        // Mise à jour de X : X = X + alpha * resid
        cblas_daxpy(*la, *alpha_rich, resid, 1, X, 1);

        // Calcul du nouveau résidu : resid = RHS - A*X
        cblas_dcopy(*la, RHS, 1, resid, 1); // resid = RHS
        cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, resid, 1);

        // Mise à jour de la norme du résidu
        norm_res = cblas_dnrm2(*la, resid, 1) / norm_rhs;

        // Sauvegarde du résidu
        resvec[*nbite] = norm_res;

        // Incrémentation du nombre d'itérations
        (*nbite)++;
    }

    // Libération de la mémoire
    free(resid);
}

void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
   for (int i = 0; i < *la; i++) {
        MB[i * (*lab)+1] = AB[i * (*lab) + 1]; 
    }
}

void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
    for (int i = 0; i < *la; i++) {
        MB[*lab * i  +1] = AB[i * (*lab) + 1]; 
        MB[*lab * i  +2] = AB[i * (*lab) + 2]; 
    }
}

void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
    
}

