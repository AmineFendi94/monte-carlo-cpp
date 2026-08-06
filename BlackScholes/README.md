# Épisode 1 à 5 — Monte Carlo sous Black–Scholes

Ce dossier accompagne la première partie de la série consacrée au pricing
Monte Carlo en C++.

## Progression

Les noms originaux des sources sont conservés. Les cinq épisodes suivent la
progression suivante :

| Épisode | Fichier associé | Contenu |
|---:|---|---|
| 1 | [`01_monte_carlo.cpp`](src/01_monte_carlo.cpp) | Présentation de la méthode Monte Carlo et application au pricing d'un call sous Black--Scholes |
| 2 | [`02_intervalle_confiance_welford.cpp`](src/02_intervalle_confiance_welford.cpp) | Calcul de l'écart-type de l'estimateur avec l'algorithme de Welford et construction de l'intervalle de confiance |
| 3 | [`03_variables_antithetiques.cpp`](src/03_variables_antithetiques.cpp) | Réduction de variance par variables antithétiques |
| 4 | [`04_variable_controle.cpp`](src/04_variable_controle.cpp) | Réduction de variance par variable de contrôle |
| 5 | [`05_antithetiques_controle_openmp.cpp`](src/05_antithetiques_controle_openmp.cpp) | Combinaison des deux méthodes de réduction de variance et parallélisation avec OpenMP |

Les épisodes 1 et 2 disposent désormais de fichiers distincts : le second
enrichit l'implémentation introduite dans le premier avec les statistiques
nécessaires à l'intervalle de confiance.

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

[Voir la série de vidéos sur YouTube](https://www.youtube.com/playlist?list=PLYDUo5ty5PuI)

## Compilation

Depuis la racine du dépôt :

```bash
cmake -S . -B build
cmake --build build
```

Les exécutables sont placés dans `build/BlackScholes/`.

Pour choisir le nombre de threads de la version OpenMP :

```bash
OMP_NUM_THREADS=4 ./build/BlackScholes/05_antithetiques_controle_openmp
```

## Support de présentation

Le fichier LaTeX utilisé pour la série est disponible dans
[`slides/MonteCarloPresentation.tex`](slides/MonteCarloPresentation.tex).
Il présente le cadre mathématique, les méthodes de réduction de variance et
la parallélisation OpenMP sans inclure de code dans les diapositives.
