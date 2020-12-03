# Ride-Sharing (Anti-)Coordination Game v1.0 (09/2020)
## Contacts
 - [David Storch](mailto:david.storch@tu-dresden.de)
 - [Malte Schröder](mailto:malte.schroeder@tu-dresden.de)

## What This Code Does
This c++-code implements a ride-sharing (anti)-coordination game as proposed in https://arxiv.org/abs/2008.11079

The present code version realizes a minimal simulation of the replicator dynamics for one set of parameters
defined in the game (S,epsilon,zeta,xi). It repeatedly solves, amongst others, a maximum weight matching problem
to pair shared ride requests. For this purpose, it embeds the 'Blossom V' algorithm proposed by Vladimier
Kolmogorov in "Mathematical Programming Computation (MPC), July 2009, 1(1): 43-67.'

## Who May Use This Code
- This code may be freely used for research purposes.
- For commercial use, please contact the authors.

## How to Run This Code
- To compile the present code, the Blossom V, v2.05, library needs to be installed.
- For research purposes, this library can be downloaded from: https://pub.ist.ac.at/~vnk/software.html
