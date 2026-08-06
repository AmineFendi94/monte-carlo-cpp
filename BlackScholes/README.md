# Épisode 1 à 5 — Monte Carlo sous Black–Scholes

Ce dossier accompagne la première partie de la série consacrée au pricing
Monte Carlo en C++.

## Progression

Les noms originaux des sources sont conservés. L'ordre conseillé pour suivre
la série est :

| Étape | Fichier | Notion principale |
|---:|---|---|
| 1 | `BlackScholes.cpp` | Simulation du call et statistiques en ligne |
| 2 | `BlackScholesAntithetic.cpp` | Variables antithétiques |
| 3 | `BlackScholesVC.cpp` | Variable de contrôle |
| 4 | `BlackScholesOpenMP.cpp` | Parallélisation OpenMP et fusion de Welford |

Le call européen sert d'exemple pédagogique : son prix analytique permet de
contrôler facilement le résultat de la simulation.

## Paramètres de référence

Les exemples utilisent les paramètres suivants :

| Paramètre | Valeur |
|---|---:|
| Prix initial $S_0$ | 197,2 |
| Strike $K$ | 200 |
| Taux sans risque $r$ | 0 |
| Volatilité $\sigma$ | 0,156582142 |
| Maturité $T$ | 1 an |

Le nombre de simulations est défini dans chaque source et peut être adapté à
la machine utilisée.

## Vidéo

Le lien vers la vidéo YouTube sera ajouté ici lors de sa publication.

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

## Support de présentation

Le fichier LaTeX utilisé pour la série est disponible dans
[`slides/MonteCarloPresentation.tex`](slides/MonteCarloPresentation.tex).
Il présente le cadre mathématique, les méthodes de réduction de variance et
la parallélisation OpenMP sans inclure de code dans les diapositives.

## Mesurer l'accélération OpenMP

Pour chaque nombre de threads $p$, mesurer le temps $T_p$, puis calculer :

$$
\operatorname{Speedup}(p)=\frac{T_1}{T_p}.
$$

| Threads | Temps $T_p$ | Speedup |
|---:|---:|---:|
| 1 | à mesurer | 1,00 |
| 2 | à mesurer | à calculer |
| 4 | à mesurer | à calculer |
| 8 | à mesurer | à calculer |
