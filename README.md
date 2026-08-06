# Série Monte Carlo en C++

Ce dépôt accompagne une série de vidéos consacrée aux méthodes de Monte Carlo
appliquées à la finance quantitative.

L'objectif est de partir d'implémentations simples, puis d'introduire
progressivement les méthodes de réduction de variance, les statistiques
numériquement stables et la parallélisation.

## Épisodes

| Dossier | Sujet |
|---|---|
| [`BlackScholes/`](BlackScholes/) | Pricing d'un call européen, variables antithétiques, variable de contrôle et OpenMP |

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

