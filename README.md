# GRN synthetic data simulation

Simulation dynamique de benchmarks de réseaux de régulation géniques pour l'entraînement d'algorithmes d'inférence.

## Table des Matières

- [Introduction](#introduction)
- [Installation](#installation)
- [Utilisation](#utilisation)
- [Tests](#tests)

## Introduction

Ce projet permet de simuler dynamiquement des réseaux de régulation géniques (GRN) en utilisant des équations différentielles ordinaires (ODEs) et de générer des benchmarks pour l'entraînement et l'évaluation d'algorithmes d'inférence. Le projet supporte trois types d'ODEs : loi Action de Masse et loi de Hill.


## Installation

### Prérequis

- Python 3.6+
- `pip`

### Étapes d'installation

1. Clonez le dépôt :
    ```sh
    git clone https://github.com/castlechild/GRNSimulation
    cd projet
    ```

2. Installez les dépendances :
    ```sh
    pip install -r requirements.txt
    ```

3. Installez le package :
    ```sh
    pip install .
    ```

## Utilisation

### Exécution de la simulation

Vous pouvez exécuter une simulation en utilisant le script principal `main.py` :

```sh
python main.py
```