# TPCalculNumérique

**Nom** : Doghmane  
**Prénom** : Amine  

---

## Compilation du TP

```bash
docker build -t poisson1d
docker run --rm -it poisson1d bash
cd app/bin
./tpPoisson1D_direct
./tp_testenv
./tpPoisson1D_iter
```

---

## Exercice 3

### Question 1

Pour utiliser BLAS et LAPACK, on déclare une matrice :


double* A;


Et on l'alloue :


A = (double*)malloc(rows * cols * sizeof(double));


Puis, pour accéder à l'élément de la matrice, on peut utiliser :


A[j * rows + i]


---

### Question 2

C’est une constante qui indique que les matrices sont stockées colonne par colonne.

---

### Question 3

La dimension principale correspond à la taille allouée pour chaque colonne ou ligne de la matrice dans la mémoire.

---

### Question 4

dgbmv effectue une opération de produit matrice-vecteur pour une matrice bande généralisée :

\[ y = \alpha A x + \beta y \]

---

### Question 5

dgbtrf effectue une factorisation LU d'une matrice bande généralisée `A` avec pivotage partiel.

---

### Question 6

dgbtrs résout un système linéaire \( A \cdot X = B \) pour une matrice bande `A`, en utilisant sa factorisation LU obtenue par `dgbtrf`.

---

### Question 7

dgbsv résout un système linéaire \( A \cdot X = B \) pour une matrice bande généralisée `A` directement en une seule étape (factorisation LU et résolution).

---

### Question 8

Comment calculer la norme du résidu relatif avec des appels BLAS ?
Rappelle répondre plus tard