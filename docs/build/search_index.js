var documenterSearchIndex = {"docs":
[{"location":"#Peel.jl-Documentation","page":"Peel.jl Documentation","title":"Peel.jl Documentation","text":"","category":"section"},{"location":"","page":"Peel.jl Documentation","title":"Peel.jl Documentation","text":"checkAndExtractSingleton","category":"page"},{"location":"#Main.PEEL.checkAndExtractSingleton","page":"Peel.jl Documentation","title":"Main.PEEL.checkAndExtractSingleton","text":"checkAndExtractSingleton(indexes,BITS;cutoff=0.000002)\n\nArguments\n\nindexes A list of arrays, being the bins and their offsets.\nBits the number of bits (being 2*qubits)\ncutoff what the entropy has to be less than\n\nChecks the 'entropy' in the bin (ie sum of the square of the fluctations around the mean of the abs value)/BITS If it is less than the cutoff admits it as a singleton. Extracts the Pauli (using majority vote if given multiple samples)\n\nNote if you wanted to implement the \"Error Correcting Code in the offsets\" as described in the appendix of the paper then this funciton would be different.\n\nReturns\n\nA tuple of: (found,index, value), where found is boolean, the index is an integer and value is float (between 0 and 1).\n\n\n\n\n\n","category":"function"},{"location":"","page":"Peel.jl Documentation","title":"Peel.jl Documentation","text":"PEEL.peelBack","category":"page"},{"location":"#Main.PEEL.peelBack","page":"Peel.jl Documentation","title":"Main.PEEL.peelBack","text":"peelBack(listOfX,listOfPs,singletonBits,singletonValue,found,ds,mappings)\n\nArguments\n\nlistOfX - the bins and their offsets <- this gets modified\nlistOfPs - the stabiliser groups used to create the bins (see generateFromPVecSamples4N)\nsingletonBits - the bit index of the Pauli to peel back\nsingletonVlaue - the value of the Pauli to peel back\nfound - a Dict of Pauli bitstrings that have been found and their value\nds - the offsets used.\nmappings - an optional mapping of qubit->qubit.\n\nPeel back algorithm tailored for noise version.\n\nCheck if we have already found the supplied Pauli (is it in found?)\n\nIf so check if the value we think we have found is greater than the previous sample\nIf not, then just return - we had already found it, and the noise has just confused us!\nIf yes, then the one we found was probably wrong and a result of the noise, just update the value.\n\n[Query this logic, it seems to work but perhaps, we may want to update the bins with the difference.]\n\nOtherwise (i.e. we hadn't found it before.)\nGet the index of that Pauli into each of the other bin sets\nAnd add or subtract (depending on relevant offset) the supplied value of the Pauli\nUpdate found by adding in this Pauli and the value.\n\n\n\n\n\n","category":"function"},{"location":"","page":"Peel.jl Documentation","title":"Peel.jl Documentation","text":"PEEL.fourPattern(n)","category":"page"},{"location":"#Main.PEEL.fourPattern-Tuple{Any}","page":"Peel.jl Documentation","title":"Main.PEEL.fourPattern","text":"fourPattern(p)\n\nArguments\n\np::Array{Array{Array{Int64,1},1},1} represents a series of commuting two qubit paulis.\n\nGiven a MUB works out what binary string will go in what bin This takes a 4 bit (two Pauli) set of commuting Paulis, eg.\n\n1-element Array{Array{Array{Int64,1},1},1}:\n [[0, 0, 0, 0], [1, 0, 1, 1], [0, 1, 1, 0], [1, 1, 0, 1]]\n\nAnd then returns the Array that shows what bin each of the 16 two qubit Paulis go into e.g.\n\n[[[\"0000\", \"0111\", \"1001\", \"1110\"], [\"0001\", \"0110\", \"1000\", \"1111\"], [\"0011\", \"0100\", \"1010\", \"1101\"], [\"0010\", \"0101\", \"1011\", \"1100\"]]]\n\nWhich tells us e.g. that \"1001\" goes into bin 0 and \"1101\" goes into bin 2.\n\nReturns\n\nan Array, containing a single array of 4 vectors each of size 4 --- splitting all 16 Paulis (rep as binary strings) into 4 bins.\n\nIf you want to think of this as a hashing function, it shows how to take 4 bits -> 2 bits using the supplied pattern.\n\n\n\n\n\n","category":"method"},{"location":"","page":"Peel.jl Documentation","title":"Peel.jl Documentation","text":"twoPattern(n)","category":"page"},{"location":"#Main.PEEL.twoPattern-Tuple{Any}","page":"Peel.jl Documentation","title":"Main.PEEL.twoPattern","text":"twoPattern(p)\n\nArguments\n\n p Array{Array{Int64,1},1} represents a series of commuting two qubit paulis.\n\nGiven a MUB works out what binary string will go in what bin This takes a 2 bit (one Pauli) set of commuting Paulis, eg.\n\n1-element Array{Array{Array{Int64,1},1},1}:\n [0, 0], [1, 1]\n\nAnd then returns a two element Array that shows what bin each of the 16 two qubit Paulis go into e.g.\n\n[\"00\", \"11\"]\n[\"01\", \"10\"]\n\nWhich tells us e.g. that \"11\" goes into bin 0 and \"01\" goes into bin 2.\n\nReturns an Array, containing an array of 2 vectors each of size 2 –- splitting all 4 Paulis (rep as binary strings) into 2 bins. If you want to think of this as a hashing function, it shows how to take 2 bits -> 1 bits using the supplied pattern.\n\n\n\n\n\n","category":"method"},{"location":"","page":"Peel.jl Documentation","title":"Peel.jl Documentation","text":"generateFromPVecSamples4N","category":"page"},{"location":"#Main.PEEL.generateFromPVecSamples4N","page":"Peel.jl Documentation","title":"Main.PEEL.generateFromPVecSamples4N","text":"generateFromPVecSamples4N(pvec,d=[],dictMapping)\n\nArguments\n\npvec The list of Paulis which are used in our experiment. Given as an array of stabilisers (see example)\nd The offset, this is the binary string 'hashed' into the Paulis given by pvec.\nmap an optional dictionary mapping qubits to other qubits.\n\nGenerates the fidelities we are going to sample from for a given offset  vector (d).\n\nBy way of example, if we supplied as a pvec\n\n[[0, 0, 0, 0], [0, 1, 1, 1], [1, 0, 0, 1], [1, 1, 1, 0]]\n\n\nThen we would get 0,7,9 and 14.\n\nIf we also suplied a d, [1 0 0 0], then we get 8, 15, 1 6\n\nA pvec as \n\n [[0, 0, 0, 0], [0, 1, 1, 1], [1, 0, 0, 1], [1, 1, 1, 0]]\n [[0, 0, 0, 0], [0, 1, 1, 1], [1, 0, 0, 1], [1, 1, 1, 0]]\n\nWould give us 0, 7, 9, 14, 112, 119, 121, 126, 144, 151, 153, 158, 224, 231, 233, 238\n\nReturns\n\nArray{Int64,1} with the relevant fidelity indexes. (zero indexed).\n\n\n\n\n\n","category":"function"},{"location":"","page":"Peel.jl Documentation","title":"Peel.jl Documentation","text":"simplexUp","category":"page"},{"location":"#Main.PEEL.simplexUp","page":"Peel.jl Documentation","title":"Main.PEEL.simplexUp","text":"simplexUp(n)\n\nArguments\n\nn (Int)  = number to project onto \n\nProjects onto the n-simplex. i.e. returns n numbers such that the sum of them is 1.\n\nReturns\n\nvector of n-numbers that add up to 1.\n\n\n\n\n\n","category":"function"},{"location":"","page":"Peel.jl Documentation","title":"Peel.jl Documentation","text":"generatePauliFromClifford","category":"page"},{"location":"#Main.PEEL.generatePauliFromClifford","page":"Peel.jl Documentation","title":"Main.PEEL.generatePauliFromClifford","text":"generatePauliFromClifford(n)\n\nArguments\n\n`cliffordPs`: Array{Float64,1}\n\nSplits the \"Clifford average\" which we recover from the protocol in Efficient Learning of quantum Noise arXiv:1907.13022. into  a fake Pauli distribution, that would average to the same Clifford average. We are going to use this to test the algorithm.\n\nWe do this by realising that each member of our 2^n distribution was in fact an average of certain Paulis. We work out how many Paulis went into the distribution (splitUp(n)) and then project the number onto an appropriate simplex (simplexUp).\n\nThis allows us to recreate a 4^n distribution that would average down to the supplied distribution. Of course there is a random  element with the projection up stage, but it is experimentally inspired if the original averaged protocol came from an experiment.\n\nReturns\n\nA Pauli probability distribution that result in the supplied probability distribution if they were the errors caused by the  avarged noise in the machine and that machine was subject to single qubit Clifford twirls.\n\n\n\n\n\n","category":"function"},{"location":"","page":"Peel.jl Documentation","title":"Peel.jl Documentation","text":"getIndexOf","category":"page"},{"location":"#Main.PEEL.getIndexOf","page":"Peel.jl Documentation","title":"Main.PEEL.getIndexOf","text":"getIndexOf(hadamardMap,no,dictMap=Dict())\n\nArguments\n\nhadamardMap: An array of hashing bins (see later)\nno: Integer  the index of the Pauli we want to know which bin it gotes into \ndictMap: Dictionary optional mapping eg qubit 2->qubit 4.\n\nGiven a map (constructed from fourPattern and twoPattern), this will  split the qubits up into the pairs/singles identified by the map and  then determine the hash destingation of a particular Pauli, i.e. which bin it goes into.  The optional dictMapping, allows us to map qubits to different qubits. This might be   relevant for different topologies. E.g. we might want the pair to span qubits 1 and 14, rather  than 1 and 2, so we would map 2 onto 14. There should be examples later.\n\nNOTE because we use the least significant (rh) qubit - there is a subte reversal of bits. See example.\n\nExample\n\nThe hadamard Map for a six qubit system might look like this\n\n4-element Array{Any,1}:\n [[\"00\", \"10\"], [\"01\", \"11\"]]\n [[\"0000\", \"0110\", \"1101\", \"1011\"], [\"1100\", \"1010\", \"0001\", \"0111\"], [\"1000\", \"1110\", \"0101\", \"0011\"], [\"0100\", \"0010\", \"1001\", \"1111\"]]\n [[\"0000\", \"1110\", \"1001\", \"0111\"], [\"1100\", \"0010\", \"0101\", \"1011\"], [\"0100\", \"1010\", \"1101\", \"0011\"], [\"1000\", \"0110\", \"0001\", \"1111\"]]\n [[\"00\", \"01\"], [\"10\", \"11\"]]\n\nWhere we have used twoPattern for qubits 0 and 5 and fourPattern for qubits 1&2 and 3&4\n\ngetIndexOf(patternAbove, 23), would then return 6.\n\nBecause 23 = \"000000010111\"  Then we map the bins as   00-> reversed to 00 -> bin 0 -> 0  0000 -> reversed to 0000 ->  bin 0 -> 00  0101 -> reversed to 1010 ->  bin 2 -> 10  11 ->  reversed to 11 -> bin 1 -> 01\n\nTherefore the index is 0 00 10 01 = 5, as we one index we return 6.\n\nReturns\n\nIndex of bin (1 indexed).\n\n\n\n\n\n","category":"function"},{"location":"","page":"Peel.jl Documentation","title":"Peel.jl Documentation","text":"constructTheKeyIndexMUB","category":"page"},{"location":"#Main.PEEL.constructTheKeyIndexMUB","page":"Peel.jl Documentation","title":"Main.PEEL.constructTheKeyIndexMUB","text":"constructTheKeyIndexMUB(indx,listOfPs,map=[])\n\nArguments\n\nindx The index of a Pauli, we want to ascertain which bins that Pauli was hashed into.\nlistOfPs The used to sample/bin the Paulis\nmap A map of qubit substitutions e.g mapping qubit 2 to 14\n\nGiven a Paulis to test, returns the index into all the bins that Paulis is hashed into  (it is assumed we have a number a sub-sampling groups). Basically just iterates over them. If we supply a map, it should be an array of maps, one for each different set of n-qubit MUBS supplied.\n\nIt uses getIndexOf to do all the heavy lifting, the documentation for that function  provides an example of the type of map it requires, this is an array of these maps.\n\nReturns and array of offsets representing the index the Pauli gets hashed to.\n\n\n\n\n\n","category":"function"}]
}
