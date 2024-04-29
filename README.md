# Timely status updates in slotted ALOHA networks with energy harvesting

This repository contains the Matlab numerical routines of the paper:

[1] K.-H. Ngo, G. Durisi, A. Munari, and F. Lazaro, and A. Graell i Amat, "Timely status updates in slotted ALOHA networks with energy harvesting," submitted to IEEE Transactions on Communications, Apr. 2024.

and of a simplified setting reported in its short version presented at Globecom 2023:

[2] K.-H. Ngo, G. Durisi, A. Graell i Amat, A. Munari, and F. Lazaro, “Age of information in slotted ALOHA with energy harvesting,” in IEEE Global Communications Conference (Globecom), Kuala Lumpur, Malaysia, Dec. 2023. URL: [https://research.chalmers.se/publication/537484/file/537484_Fulltext.pdf](https://research.chalmers.se/publication/537484/file/537484_Fulltext.pdf).

Please cite the aforementioned papers if you use this code.

## Content of the repository

This repository contains the following main files used for [1]:

1. `exact.m`: evaluate the exact average age of information (AoI) via a first-order Markov analysis.
2. `simulation.m`: evaluate the average AoI, AVP, and throughput via a simulation of the protocol.
3. `approximation.m`: evaluate approximations of the average AoI, AVP, and throughput obtained via the phase-type distribution.
4. `optimize_ptx.m`: optimize the transmission probabilities to minimize the average AoI, minimize the AVP, or maximize the throughput.

The folder `simplified_setting` contains the corresponding files used for the setting in [2].

The folder `helpers` contains the following supplementary functions:
1. `compute_w_capture.m`: evaluate the successful decoding probability for each device when the decoding is with capture.
2. `interger_partitions.m`: partition a given interger into multiple integers.
3. `uniqueperms.m`: generate all unique permutations of elements in a vector.
4. `randi_distr.m`: generate a matrix with random entries distributed according to a specified discrete distribution.
5. `fminsearchcon.m`: optimization using the Nelder-Mead method; extension of FMINSEARCHBND in Matlab with general inequality constraints.

The functions `interger_partitions`, `uniqueperms`, and `fminsearchcon` were written by other people. The references are provided in the files.
