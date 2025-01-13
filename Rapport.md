# Rapport de TD/TP : Calcul Numérique

## Introduction

L’équation de la chaleur est un problème central en physique et en mathématiques appliquées, décrivant la diffusion de la température dans un milieu homogène immobile. Ce TD/TP se concentre sur la résolution numérique de l’équation stationnaire en 1D sans terme source (« g = 0 ») par discrétisation avec un schéma de différences finies centrées d’ordre 2.

La solution analytique est donnée par :
\[
T(x) = T_0 + x(T_1 - T_0),
\]
ou \(T_0\) et \(T_1\) sont les températures aux bords.

Ce rapport présente les étapes de résolution à travers des méthodes directes et itératives, en utilisant des bibliothèques de calcul BLAS et LAPACK. Le travail s’est déroulé en trois parties : mise en place de l’environnement, méthodes directes et méthodes itératives.

---

## Partie 1 : Préparation de l'environnement

### Mise en place

L'environnement de travail est configuré via un conteneur Docker permettant l'utilisation des bibliothèques CBLAS et LAPACK. Une série de tests a permis de valider leur bonne installation. Voici les résultats principaux :

- Valeur exponentielle (« e ») : 2.718282.
- Valeurs maximales en simple et double précision :
  - flt\_max = 3.402823e+38
  - dbl\_max = 1.797693e+308
- Précision relative :
  - flt\_epsilon = 1.192093e-07
  - dbl\_epsilon = 2.220446e-16

Les fonctions BLAS comme `DCOPY` ont été testées avec succès, confirmant la présence d'un environnement fonctionnel.


## Instructions de Compilation et d’Exécution

Pour compiler et exécuter les programmes, suivez les étapes suivantes :

1. Construisez l’image Docker :
   ```bash
   docker build -t poisson1d .
   ```

2. Lancez un conteneur Docker interactif :
   ```bash
   docker run --rm -it poisson1d bash
   ```

3. Accédez au répertoire des exécutables :
   ```bash
   cd app/bin
   ```

4. Exécutez les différents programmes :
   - **Test de l’environnement :**
     ```bash
     ./tp_testenv
     ```
   - **Méthodes directes :**
     ```bash
     ./tpPoisson1D_direct
     ```
   - **Méthodes itératives :**
     ```bash
     ./tpPoisson1D_iter
     ```

---

## Partie 2 : Méthodes Directes



### Structure

1. **Méthode directe et stockage bande**
   - Préparation des données et construction de la matrice
   - Utilisation des fonctions LAPACK
   - Résultats obtenus
2. **Validation et comparaison des méthodes**
   - Validation du stockage GB avec `dgbmv`
   - Comparaison avec les autres approches
3. **Validation de l’environnement et tests des exécutables**

### Méthode Directe et Stockage Bande

#### Préparation des données et construction de la matrice

La matrice associée au problème est tridiagonale et représentée dans un format « General Band » (GB) pour optimiser la mémoire et les calculs. Cette matrice est construite à l’aide de la fonction `set_GB_operator_colMajor_poisson1D`. Elle organise les éléments diagonaux et hors-diagonaux dans un tableau compact pour éviter de stocker les zéros inutiles.




Le vecteur des seconds membres (« RHS ») est généré par la fonction `set_dense_RHS_DBC_1D`. Ce vecteur incorpore les conditions aux limites du problème : les températures imposées \(T_0\) et \(T_1\) sont intégrées dans les premières et dernières entrées respectivement, tandis que les autres valeurs sont initialisées à zéro.

Le maillage 1D est également construit grâce à `set_grid_points_1D`, qui génère un ensemble de points uniformément espacés dans le domaine considéré. Ces points servent de base pour la discrétisation spatiale et permettent de calculer la solution analytique.


![](images/grid_points.png)


Détails supplémentaires sur le stockage bande : Le stockage GB réduit significativement l’occupation mémoire en ne sauvegardant que les diagonales pertinentes. Pour une matrice tridiagonale de taille \(n\), seule une bande de largeur \(3\) (diagonale principale et deux diagonales adjacentes) est conservée, au lieu de la matrice complète \(n \times n\). Cela permet d’économiser de la mémoire et d’accélérer les opérations.

#### Utilisation des fonctions LAPACK

- **Factorisation LU :** La routine `dgbtrf` décompose la matrice GB en deux matrices triangulaires \(L\) (inférieure) et \(U\) (supérieure). Cette factorisation est essentielle pour résoudre efficacement le système linéaire associé.

- **Résolution :** Une fois la matrice factorisée, `dgbtrs` est utilisée pour obtenir la solution numérique. Elle exploite les matrices \(L\) et \(U\) pour résoudre rapidement le système sans recalculer la factorisation.

- **Validation :** Les solutions numériques sont comparées à la solution analytique (« EX_SOL.dat ») à l’aide de l’erreur relative. Cette métrique évalue la précision de la méthode en mesurant l’écart entre la solution exacte et la solution calculée.

#### Résultats obtenus

- **Précision :** Les solutions numériques sont proches de la solution analytique, avec une erreur relative typiquement de l’ordre de \(10^{-16}\). Cela montre que les méthodes directes sont capables de produire des résultats très précis pour les matrices de taille modérée.

- **Validations :** Les fichiers générés, tels que `LU.dat`, confirment que la factorisation LU est correcte. Les matrices \(L\) et \(U\) respectent la structure tridiagonale initiale.



- **Performances :** Les temps d’exécution pour `dgbtrf` et `dgbtrs` ont été mesurés. La factorisation présente une complexité temporelle en \(O(n^3)\), tandis que la résolution est plus rapide, en \(O(n^2)\).

![solution_comparison](images/solution_comparison.png)


Un graphique comparatif met en évidence la concordance entre les solutions analytiques et numériques. Ce graphique permet de visualiser directement l’efficacité de la méthode.

### Validation et Comparaison des Méthodes

#### Validation du Stockage GB avec `dgbmv`

- **Produit matrice-vecteur :** La fonction `dgbmv` est utilisée pour effectuer un produit matrice-vecteur avec la matrice GB. Ce test valide que les éléments de la matrice sont correctement stockés et accessibles dans leur format compact.

- **Résultats :** Les tests montrent que le produit est exact, avec un résidu proche de zéro. Cela confirme que la matrice a été correctement construite et que les calculs sont fiables.

#### Comparaison avec d’autres approches

- **Avantages des méthodes directes :** Les méthodes directes, comme la factorisation LU, offrent une précision élevée et une stabilité numérique. Elles conviennent particulièrement pour les matrices de taille modérée.

- **Limites :** Leur principal inconvénient est leur complexité temporelle et spatiale, qui augmente rapidement avec la taille de la matrice. Cela les rend moins adaptées pour les grands systèmes.

- **Perspectives :** Les méthodes itératives, telles que Richardson ou Gauss-Seidel, représentent une alternative intéressante pour les grands systèmes. Ces méthodes exploitent la structure creuse des matrices pour réduire les coûts de calcul.

### Validation de l’Environnement et Tests des Exécutables

Avant d’exécuter les méthodes de résolution, l’environnement de calcul a été testé pour valider la configuration des bibliothèques et s’assurer de la précision des calculs.

#### Résultats de Test de l’Environnement

1. **Précision des Calculs :**
   - Les valeurs importantes comme l’exponentielle (\(e\)), les limites de précision en simple et double précision, ainsi que les epsilon machine, ont été vérifiées. Ces résultats confirment que l’environnement est correctement configuré pour des calculs numériques fiables.
   
   ```
   The exponantial value is e = 2.718282
   The maximum single precision value from values.h is maxfloat = 3.402823e+38
   The epsilon in double precision value is dbl_epsilon = 2.220446e-16
   ```

2. **Tests BLAS/LAPACK :**
   - Une copie vecteur à vecteur (`DCOPY`) et des tests de base sur des vecteurs montrent que les fonctions BLAS fonctionnent correctement.

   ```
   x[0] = 1.000000, y[0] = 6.000000
   Test DCOPY y <- x
   y[0] = 1.000000
   ```

3. **Exécutions des Solutions :**
   - Pour le cas direct, l’erreur relative est extrêmement faible (\(2.393518e-16\)), validant la précision des solutions calculées avec LAPACK.

   ```
   The relative forward error is relres = 2.393518e-16
   ```

- Lors d’autres configurations, des erreurs relatives plus élevées (\(6.813851e-01\)) indiquent des cas mal conditionnés ou des solutions moins précises, ce qui nécessite une analyse plus poussée.

### Annexes

1. **Fonctions clés utilisées :**
   - `set_GB_operator_colMajor_poisson1D` : Construction de la matrice au format GB.
   - `set_dense_RHS_DBC_1D` : Génération du vecteur RHS avec conditions aux limites.
   - `set_grid_points_1D` : Construction du maillage 1D pour la discrétisation spatiale.

2. **Glossaire des fichiers :**
   - `RHS.dat` : Contient le vecteur des seconds membres.
   - `EX_SOL.dat` : Contient la solution analytique.
   - `X_grid.dat` : Contient les coordonnées des points du maillage.
   - `LU.dat` : Contient les matrices \(L\) et \(U\) après factorisation.

Ces annexes permettent de clarifier les rôles des différentes fonctions et fichiers dans le cadre de ce TP.

---

## Partie 3 : Méthodes Itératives

### Algorithme de Richardson

L'algorithme de Richardson est implémenté avec un paramètre \( \alpha \) optimisé, calculé à partir des valeurs propres de la matrice. La convergence est analysée par le suivi de la norme du résidu.


![convergence richardson](images/convergence.png)


### Résultats

- L'erreur relative par rapport à la solution analytique est de \( 5.093853 \times 10^{-3} \).
- La convergence est rapide pour les matrices tridiagonales.

---

## Analyse Comparative

| Critère            | Méthode Directe           | Méthodes Itératives       |
|--------------------|---------------------------|---------------------------|
| **Précision**      | Très élevée (écart faible)  | Modérée                  |
| **Complexité**     | \( O(n^3) \) pour `dgbtrf`  | Linéaire (à chaque itération) |
| **Utilisation**    | Systèmes de taille modérée | Grandes matrices          |

---

## Conclusion

Ce rapport met en avant deux approches complémentaires pour la résolution de l’équation de Poisson en 1D :

1. **Méthodes Directes :** Idéales pour les problèmes de taille modérée, avec une précision très élevée.
2. **Méthodes Itératives :** Plus adaptées aux grands syst
