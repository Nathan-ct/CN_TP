# Rapport TP - Calcul Numérique






## Exercice 3. Référence et utilisation de BLAS/LAPACK

À l’aide de la documentation de BLAS et LAPACK, répondez aux questions suivantes :

1. **En C, comment doit-on déclarer et allouer une matrice pour utiliser BLAS et LAPACK ?**  
   Pour utiliser BLAS et LAPACK en C, il faut déclarer ses matrices comme des tableaux 1D, les allouer dynamiquement, et gérer l'accès avec le format column-major.

2. **Quelle est la signification de la constante `LAPACK_COL_MAJOR` ?**  
   La constante `LAPACK_COL_MAJOR` indique que la matrice est stockée en format column-major.

3. **À quoi correspond la dimension principale (leading dimension), généralement notée `ld` ?**  
   La dimension principale correspond à la taille mémoire allouée pour une colonne de la matrice.

4. **Que fait la fonction `dgbmv` ? Quelle méthode implémente-t-elle ?**  
   La fonction `dgbmv` réalise une multiplication matrice-vecteur avec une matrice bande.

5. **Que fait la fonction `dgbtrf` ? Quelle méthode implémente-t-elle ?**  
   La fonction `dgbtrf` est une fonction LAPACK qui effectue la décomposition LU d'une matrice bande avec pivotage partiel. Elle utilise une version adaptée de l'algorithme LU pour les matrices bandes afin de minimiser les calculs et la mémoire.

6. **Que fait la fonction `dgbtrs` ? Quelle méthode implémente-t-elle ?**  
   La fonction `dgbtrs` est une fonction LAPACK qui résout un système linéaire de la forme : \( AX = B \).  
   Elle utilise la décomposition LU calculée par `dgbtrf` pour résoudre le système linéaire en effectuant une substitution directe avant/arrière.

7. **Que fait la fonction `dgbsv` ? Quelle méthode implémente-t-elle ?**  
   La fonction `dgbsv` résout un système linéaire avec une matrice bande.  
   Elle effectue la décomposition LU (`dgbtrf`) suivie de la substitution avant/arrière (`dgbtrs`) pour résoudre le système. Elle combine donc les deux étapes.

8. **Comment calculer la norme du résidu relatif avec des appels BLAS ?**  
   La norme du résidu relatif est obtenue en utilisant `dgemv` et `dnrm2`.

---

## Exercice 5. `DGBTRF`, `DGBTRS`, `DGBSV`

Pour une matrice tridiagonale de taille \( n \times n \) :

- **DGBTRF (Factorisation LU bande)** :
  - **Complexité** : \( O(n) \) opérations car la matrice est tridiagonale.
  - **Stockage** : \( O(3n) \), car nous avons 3 diagonales.

- **DGBTRS (Résolution triangulaire)** :
  - **Complexité** : \( O(n) \) opérations.
  - Utilise les matrices \( L \) et \( U \) stockées dans \( AB \).

- **DGBSV (Combinaison des deux)** :
  - **Complexité totale** : \( O(n) \) opérations.
  - Plus efficace que d'appeler `DGBTRF` et `DGBTRS` séparément car optimisé.

### Les avantages de chaque méthode :
- **DGBTRF + DGBTRS** : permet de réutiliser la factorisation pour plusieurs seconds membres.
- **DGBSV** : plus simple à utiliser et légèrement plus efficace pour un seul second membre.


## Exercice 7-8-9. Implémentation C

### Principe des méthodes

### a) Méthode de Richardson
- Formule : `x(k+1) = x(k) + α(b - Ax(k))`
- α est un paramètre de relaxation
- La convergence dépend du choix de α

### b) Méthode de Jacobi
- Utilise la décomposition `A = D + R`
- Formule : `x(k+1) = D⁻¹(b - Rx(k))`
- Pour Poisson 1D : équivalent à Richardson avec α = 0.5

### c) Méthode de Gauss-Seidel
- Utilise les valeurs mises à jour immédiatement
- Plus rapide que Jacobi car elle utilise l'information plus récente

### Implémentation

Les trois méthodes ont été implémentées dans :
- `lib_poisson1D_richardson.c` pour Richardson
- Ajout de `jacobi_poisson1D()` pour Jacobi
- Ajout de `gauss_seidel_poisson1D()` pour Gauss-Seidel

### Tests et Résultats

#### Configuration des tests

- Tolérance : `1.000000e-03`
- Nombre maximum d'itérations : variable
- Matrice tridiagonale (-1, 2, -1)

#### Résultats comparatifs

#### a) Méthode de Richardson (IMPLEM=0)
- Alpha optimal : `0.500000`
- Convergence en `126` itérations
- Résidu initial : `1.000000e+00`
- Évolution des premiers résidus :
  - Iteration 1 : `5.000000e-01`
  - Iteration 2 : `3.535534e-01`
  - Iteration 3 : `2.795085e-01`
  - Iteration 4 : `2.338536e-01`
  - Iteration 5 : `2.025231e-01`
- Résidu final : `9.669311e-04`

#### b) Méthode de Jacobi (IMPLEM=1)
- Alpha optimal : `0.500000`
- Convergence : `125` itérations
- Résidu initial : `5.000000e-01`
- Évolution des premiers résidus :
  - Iteration 1 : `3.535534e-01`
  - Iteration 2 : `2.795085e-01`
  - Iteration 3 : `2.338536e-01`
  - Iteration 4 : `2.025231e-01`
  - Iteration 5 : `1.795176e-01`
- Résidu final : `9.669311e-04`

#### c) Méthode de Gauss-Seidel (IMPLEM=2)
- Alpha optimal : `0.500000`
- Convergence : `63` itérations
- Résidu initial : `4.893320e-01`
- Évolution des premiers résidus :
  - Iteration 1 : `2.751859e-01`
  - Iteration 2 : `1.874051e-01`
  - Iteration 3 : `1.457831e-01`
  - Iteration 4 : `1.222933e-01`
  - Iteration 5 : `1.069563e-01`
- Résidu final : `9.770070e-04`

###  Analyse et Conclusions

Les résultats obtenus montrent les différences en termes de convergence et d'efficacité entre les trois méthodes implémentées pour résoudre l'équation de Poisson 1D.

### Méthode de Richardson
- Bien que cette méthode soit simple à implémenter, sa convergence est relativement lente avec `126` itérations nécessaires pour atteindre la tolérance spécifiée.
- L'évolution du résidu montre une décroissance régulière mais plus lente par rapport aux autres méthodes.
- Le choix optimal du paramètre α=`0.5` est crucial pour assurer la convergence.

### Méthode de Jacobi
- Cette méthode, mathématiquement équivalente à Richardson pour α=`0.5`, converge en `125` itérations, confirmant son comportement similaire dans ce cas spécifique.
- La structure itérative basée sur la décomposition de la matrice est utile pour des implémentations parallèles, mais ici elle n’apporte pas de gain en convergence par rapport à Richardson.

### Méthode de Gauss-Seidel
- Cette méthode converge significativement plus vite que les deux autres, avec seulement `63` itérations nécessaires.
- L’utilisation des valeurs mises à jour à chaque étape améliore la vitesse de convergence, ce qui en fait une méthode particulièrement efficace pour les matrices bien conditionnées comme dans le cas de Poisson 1D.

### Comparaison des performances

- La méthode de Gauss-Seidel est clairement la plus performante parmi les trois, divisant presque par deux le nombre d'itérations nécessaires par rapport à Jacobi et Richardson.
- La méthode de Richardson est la moins efficace, mais elle reste intéressante en tant qu’approche de base et pour étudier l’effet du paramètre de relaxation α.
- La méthode de Jacobi, bien que plus performante que Richardson dans des scénarios plus complexes (notamment en parallélisation), ne présente pas d’avantage significatif dans ce cas précis.



