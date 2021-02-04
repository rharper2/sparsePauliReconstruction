## Introduction

This repository serves to hold the code for the paper **Fast estimation of sparse quantum noise** by *Harper, Yu and Flammia* [ArXiv:2007.07901](https://arxiv.org/abs/2007.07901).

Full documentation can be found at [https://rharper2.github.io/sparsePauliReconstruction/docs/build/index.html](https://rharper2.github.io/sparsePauliReconstruction/docs/build/index.html).

## Summary

This series of files is not just the code used in the paper, but it also provides walkthoughs (in the way of Julia workbooks) for many of the core concepts used in the paper in the hope that his might be of relevance and interest to those exploring the world of characterizing the noise in quantum devices and verifying the operation of such devices.

All the work books are saved as jupyter notebooks, with html versions. Although github does a great job of rendering IPython notebooks, some of the software used prints arrays in latex form (so they look good in the workbook) and github is struggling with that (it shows the raw latex). So for online browsing I would recommend using the html version. Better still download the repo, fire up your own jupyter notebook server and view the workbooks (interactively).

## Scalable Estimation

This is the main workbook which uses data taken from the IBM Quantum Experience (Melbourne device, when it only had 14 qubits), uses that to create a full Pauli distribution and then attempts to reconstruct the distribution from limited sampling of the eigenvalues corrupted by varying levels of noise. It contains the code, analysis and figures that appear in the paper. It does, however, assume a certain level of knowledge which is the point of the workbooks below.

It implicilty uses the algorithm detailed in **Efficient Learning of quantum Noise** [arXiv:1907.13022](https://arxiv.org/abs/1907.13022) , code for which is located at https://github.com/rharper2/Juqst.jl. Python code to run such experiments on the IBM Quantum Experience (using qiskit) can be found on https://github.com/rharper2/query_ibmq.

## Hadamard Basics and Observations

This workbook is a basic (gentle) introduction to the idea of Pauli channels/Pauli erorr probabilities/Pauli eigenvalues and how the Hadamard transform converts between them. It goes through a lot of the notation we use in the paper in detail. Many will be familiar with what is contained in the workbook - but for some it may help provide the key to understanding later workbooks.

## Scalable Estimation - Basic Concepts

This workbook  uses what we learnt above (there is some repetition at the beginning) to show exactly how the PEELing decoder works in a toy 2 qubit system. It sets up a fake Pauli distribution, assumes we can sample without noise and shows how by sampling the eigenvalues we can operate the decoder to calcluate the global Pauli error probabilities. Obviously the real algorithm, which has to deal with noise, is slightly more complicated and clearly a two qubit system is almost tediously trivial - but the idea is to use a system small enough to easily manipulate to introduce the major concepts behind the algorithm.

## Scalable Estimation - Experimental Basics

This workbook goes through the basic ideas behind experimental design, for the experiments desgined to extract the eigenvalues we need. Here we use a trivial 6 qubit system. 

### Initial things to note about the experimental example

This is of marginal utility of only a 6 qubit system! There are only 4096 Paulis to measure in a 6 qubit system. The protocol requires a minimum number of 4n+2 experiments, each measuring ($2^6$) 64 possible outcomes. 26*64 = 1664 measurements - we can reconstruct all the Paulis for not much more than this! For this type of protocol we are assuming a $\delta$ of about $0.25$, ie  $4^{0.25\times6}$, so we are expecting that we will only be seeking to recover the 8-9 highest weight Paulis - obviously as system size increases that is where the algorithm shines. So for instance with the 14 qubit example, in a later workbook - it makes a lot more sense.

---
<div style="border: 3px solid red"><p style="padding:5px 10px 5px 10px;"><strong>Warning:</strong> As mentioned above the "Experimental Basics" workbook is a long and tedious workbook that is, ulitmately, in some sense disappointing. Because we simulate the entire experiment we can (and do) look at all the intermediate results. Six qubits takes about a minute to simulate a run, and we need to do 4*6+2 of them. At then end we get 8-14 numbers that are *approximately* correct. I'll try and include a qiskit simulation, which will have a less realistic noise model - but will be more qubits. The real use of the algorithm is shown in the work book which was used for the paper (the Scalabable Estimation workbook). The Scalable Estimation workbook is based off actual 14 qubit experimental data and recovers a lot more information (because it is a bigger system). But if you want the whole gory detail, the system here is small enough to grasp everything that is happening but big enough to be not entirely trivial, although it's a lot of work for very little (because of the small size of the system).</p></div>

---

