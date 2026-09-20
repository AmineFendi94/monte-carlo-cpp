# Série Monte Carlo en C++

<p align="center">
  <img src="assets/cover.png" alt="Monte Carlo en C++ — pricing d'options" width="100%">
</p>

<p align="center">
  <a href="https://www.youtube.com/playlist?list=PLYDUo5ty5PuI"><strong>▶ Voir la série complète sur YouTube</strong></a>
</p>

Ce dépôt accompagne une série de vidéos consacrée aux méthodes de Monte Carlo
appliquées à la finance quantitative.

L'objectif est de partir d'implémentations simples, puis d'introduire
progressivement les méthodes de réduction de variance, les statistiques
numériquement stables et la parallélisation.

## Épisodes

| Dossier | Sujet |
|---|---|
| [`BlackScholes/`](BlackScholes/) | Pricing de calls européens et barrières, réduction de variance et OpenMP |

## Construire les exemples

Prérequis :

- un compilateur compatible C++17 ;
- CMake 3.16 ou plus récent ;
- OpenMP.

```bash
cmake -S . -B build
cmake --build build
```

## Organisation

Chaque dossier correspond à un thème ou à un épisode. Il contient ses sources,
sa documentation et, lorsque cela est utile, les supports associés.

## Licence

Le code et les supports de ce dépôt sont distribués sous
[licence MIT](LICENSE).
