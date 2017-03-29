Full Configuration Interaction
==============

Mini Project 1<br />
MPhil in Scientific Computing<br />
Churchill College<br />
University of Cambridge

*This work was produced as part of the miniproject requirements of the MPhil in Scientific Computing course I undertook as a student of Churchill College, University of Cambridge. This work was done under the supervision of Dr. Alex Thom and with funding from the Sir Winston Churchill Foundation of the USA.*

## Introduction

This program is designed to calculate the FCI energies of a system using precalculated integrals. The Hamiltonian is set up and stored in a basis of Slater determinants and then diagonalized using the Davidson diagonalization algorithm. The program naturally interfaces with the [Q-CHEM](http://www.q-chem.com/) electronic structure package, utilizing the integrals outputted from Q-CHEM.

## Compiling

The FCI program uses the [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) library for linear algebra. The program also utilizes the C++11 standard and has parallelization implemented using [OpenMP](http://www.openmp.org/). Note that OpenMP cannot be disabled as timing is done through OpenMP functions. The following should allow for compilation of the program.

```
% g++ -std=c++11 -I /path/to/eigen Hamiltonian.cpp Diagonalization.cpp ReadInput.cpp -O3 -fopenmp -o FCI
```

## Running the Program

The program takes two command line inputs. These are, in order, the input filename and the output filename. The first command line input is the file that contains all the settings and values of the integrals. The second command line input is the filename of the output. Alternatively the program can be simply ran without any command line inputs. A prompt will ask for these to be input individually.

## Format of the Input File

The input file is formatted as such. First, a string of input parameters should be listed. In order, they are:
- (Integer, Positive) The number of alpha electrons.
- (Integer, Positive) The number of alpha spin orbitals.
- (Integer, Positive) The number of beta electrons.
- (Integer, Positive) The number of beta spin orbitals.
- (Integer, Positive) The number of FCI solutions desired, ordered from lowest to highest.

Next, the values for the two electron integrals (nm|kl) are listed in the following format
```
(nm|kl)     n     m     k     l
```
It should be noted that the nuclear repulsion has n, m, k, and l set to zero and the one electron integrals are labelled by n and m while k and l are set to zero. This is the format of Q-Chem. An example for H<sub>2</sub> is shown below. In this example, there is 1 alpha electron and 1 beta electron. The space is spanned by 4 alpha spin orbitals and 4 beta spin orbitals. The two lowest solutions are desired.
```
1 4 1 4 2
  0.64985185942031   1   1   1   1
  0.16712550470738   1   3   1   1
  0.080102886434995  1   2   1   2
  0.07936780580498   1   4   1   2
  (And the rest of the integrals)
```
