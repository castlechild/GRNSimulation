#!/usr/bin/env python3
import numpy as np


def stochastiqueNoise(production, degradation, noise_amplitude):
    return noise_amplitude * (
        np.sqrt(np.abs(production)) * np.random.normal()
        + np.sqrt(np.abs(degradation)) * np.random.normal()
        )


def main():
    pass


if __name__ == "__main__":
    main()
