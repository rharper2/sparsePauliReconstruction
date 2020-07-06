# Scalable Estimation Documentation



```@contents
Pages = ["peel.md", "localPeel.md"]
Depth = 2
```

# Introduction 

This is the documentation of the code used to implement the algorithm discussed in the paper **Fast estimation of sparse quantum noise** by *Harper, Yu and Flammia* (in production).

There are a number of IJulia workbooks that accompany this code that detail the use of the software and the implementation of the algorithm.

They are as follows:
- one
- two 
- three
- fourt


## Scalable Estimation

The main workbook uses data taken from the IBM Quantum Experience (Melbourne device, when it only had 14 qubits), uses that to create a full Pauli distribution and then attempts to reconstruct the distribution from limited sampling of the eigenvalues corrupted by varying levels of noise. It contains the code, analysis and figures that appear in the paper. It does, however, assume a certain level of knowledge which is the point of the workbooks also contained in this repository

It implicilty uses the algorithm detailed in **Efficient Learning of quantum Noise** [arXiv:1907.13022](https://arxiv.org/abs/1907.13022) , code for which is located at https://github.com/rharper2/Juqst.jl. Python code to run such experiments on the IBM Quantum Experience (using qiskit) can be found on https://github.com/rharper2/query_ibmq.


Figure 2 form the paper, shows the type of recovery that is possible using the code and cicuits discussed here

![Figure](./figure2.png)


Copyright: Robin Harper 2019-2020

