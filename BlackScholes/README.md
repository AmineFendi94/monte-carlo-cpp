# Épisode 1 — Monte Carlo et Black–Scholes

Ce dossier accompagne la première partie de la série consacrée au pricing
Monte Carlo en C++.

## Progression

Les fichiers de `src/` présentent progressivement :

1. `BlackScholes.cpp` : simulation du call européen et statistiques en ligne ;
2. `BlackScholesAntithetic.cpp` : ajout des variables antithétiques ;
3. `BlackScholesVC.cpp` : ajout d'une variable de contrôle ;
4. `BlackScholesOpenMP.cpp` : combinaison des méthodes et parallélisation OpenMP.

Le call européen sert d'exemple pédagogique : son prix analytique permet de
contrôler facilement le résultat de la simulation.

## Compilation

Depuis la racine du dépôt :

```bash
cmake -S . -B build
cmake --build build
```

Les exécutables sont placés dans `build/BlackScholes/`.

Pour choisir le nombre de threads de la version OpenMP :

```bash
OMP_NUM_THREADS=4 ./build/BlackScholes/black_scholes_openmp
```

