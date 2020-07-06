using Juqst

potentialSingles = [
                    [[0,0],[0,1]], # IX
                    [[0,0],[1,0]], # IY
                    [[0,0],[1,1]], # IZ
                    ]


all2QlMuBs =  [  [[0,0,0,0],[1,1,0,1],[1,0,1,1],[0,1,1,0]], #II ZX YZ XY
                 [[0,0,0,0],[1,1,1,0],[0,1,1,1],[1,0,0,1]], #II ZY XZ YX
                 [[0,0,0,0],[0,0,0,1],[0,1,0,0],[0,1,0,1]], #II IX XI XX
                 [[0,0,0,0],[0,0,1,0],[1,0,0,0],[1,0,1,0]], #II IY YI YY
                 [[0,0,0,0],[0,0,1,1],[1,1,0,0],[1,1,1,1]]] #II IZ ZI ZZ

# We only want to choose the first two types for the initial runs
potentialMuBs = [all2QlMuBs[1],all2QlMuBs[2]]


pX = [0 1;1 0]
pY = [0 -im;im 0]
pZ = [1 0;0 -1]
pI = [1 0;0 1]
paulis = [pI,pX,pY,pZ]
superPaulis = [makeSuper(x) for x in paulis];

"""
generateSensibleSubsamples(experiments)
## Arguments
	experiments: A list of tuples, showing how to locally divide the qubits.
Each experiment needs to be for the same number of qubits.
## Example
If you pass in [(2,2,2),(1,2,2,1)]
Then this would generage the needed data for two expeiments, each for 6 qubits.
- The first would use two qubit gates on qubits 0&1, 2&3 and 4&5.
- The second would use single qubit gates on qubits 0 and 5, and two qubit gates on 1&2 and 3&4.

With three qubits, we might call it thus:
```
(experiments,paulisAll,ds) = generateSensibleSubsamples([(2,1),(1,2)])
print("$experiments\n")
print("$paulisAll\n")
```

```
Any[Any[(2, 1), (1, 1)], Any[(1, 3), (2, 2)]]
Any[Any[[[0, 0, 0, 0], [1, 1, 0, 1], [1, 0, 1, 1], [0, 1, 1, 0]], [[0, 0], [0, 1]]], Any[[[0, 0], [1, 1]], [[0, 0, 0, 0], [1, 1, 1, 0], [0, 1, 1, 1], [1, 0, 0, 1]]]]
```

The Experiment numbers are defined like this:

```
potentialSingles = [
                    [[0,0],[0,1]], # IX
                    [[0,0],[1,0]], # IY
                    [[0,0],[1,1]], # IZ
                    ]


all2QlMuBs =  [  [[0,0,0,0],[1,1,0,1],[1,0,1,1],[0,1,1,0]], #II ZX YZ XY
                 [[0,0,0,0],[1,1,1,0],[0,1,1,1],[1,0,0,1]], #II ZY XZ YX
                 [[0,0,0,0],[0,0,0,1],[0,1,0,0],[0,1,0,1]], #II IX XI XX
                 [[0,0,0,0],[0,0,1,0],[1,0,0,0],[1,0,1,0]], #II IY YI YY
                 [[0,0,0,0],[0,0,1,1],[1,1,0,0],[1,1,1,1]]] #II IZ ZI ZZ
```


## Returns
The following tuple (chosenExperiments, chosenPaulis, ds), where
- chosenExperiments is a list of experiment list, where each experiment list contains a tuple of qubit size and experiment number.
- chosenPaulis are the Paulis that will be generated by each tuple in an experiment list (using the 01->X 10->Y 11->Z mapping) 
- ds is just a convenient bit-array of offsets of 2*noOfQubits, used in working out offsets for the decoder.
"""
function generateSensibleSubsamples(experiments)
	qubits = sum(experiments[1])
    for i in experiments
        @assert(qubits == sum(i),"All experiments have to use the same number of qubits")
    end
    chosenExperiments = []
    chosenPaulis = []
    for e in experiments
        thisExperiment =[]
        thisPaulis = []
        for q in e
            if q ==2
                chosen = rand(1:2)
                push!(thisExperiment,(2,chosen))
                push!(thisPaulis,all2QlMuBs[chosen])
            elseif q == 1
                chosen = rand(1:3)
                push!(thisExperiment,(1,chosen))
                push!(thisPaulis,potentialSingles[chosen])
            else
                @assert(false,"We only deal with one and two qubit splis just now\n")
            end
        end
        push!(chosenExperiments,thisExperiment)
        push!(chosenPaulis,thisPaulis)
    end
    ds = vcat(
        [[0 for _ = 1:2*qubits]],
        [map(x->parse(Int,x),collect(reverse(string(i,base=2,pad=2*qubits)))) for i in [2^b for b=0:2*qubits-1]]);
    return (chosenExperiments,chosenPaulis,ds)

end

"""
addCircuit(circuit,experiment_type,qubits_id)

## Arguments
   circuit - the qiskit circuit we will add gates to.
   experiment-type: (Integer,Integer) - the stabilizer group we wont to generate (see below)
   qubits_id: Array of Int - the 'number' of the qubits the circuit is to be on.

The stabilizer circuits vary depending on the experiment-type. This is a tuple, the first number indicating the number of 
qubits (currently 1 or 2), the second which experiment currently 1:3 for single qubits and 1:5 for two qubits.

The qubits_id will be used to work out which qubits in the circuit to use. They are "JULIA INDEXED".

The qubits are identified via the circuit ```qubits = circuit.qubits```. But because of the way PyCall works, they will be 
in a 1 indexed array. e.g. qiskitits qubit_0 will be qubits[1].

Therefore if you want the circuit to be on qubits 0 and 1, you pass in [1,2] into qubits_id.

The circuits that will be added will generate the following for two qubits its:

1.    [[0,0,0,0],[1,1,0,1],[1,0,1,1],[0,1,1,0]], #II ZX YZ XY
2.    [[0,0,0,0],[1,1,1,0],[0,1,1,1],[1,0,0,1]], #II ZY XZ YX
3.    [[0,0,0,0],[0,0,0,1],[0,1,0,0],[0,1,0,1]], #II IX XI XX
4.    [[0,0,0,0],[0,0,1,0],[1,0,0,0],[1,0,1,0]], #II IY YI YY
5.    [[0,0,0,0],[0,0,1,1],[1,1,0,0],[1,1,1,1]]] #II IZ ZI ZZ

For one qubit its

1.    [[0,0],[0,1]], # IX
2.    [[0,0],[1,0]], # IY
3.    [[0,0],[1,1]], # IZ

## Returns nothing

The qiskit circuit will have had the gates added.
"""
function addCircuit(circuit,experiment_type,qubits_id)
    qubits = circuit.qubits
    (q,no) = experiment_type
    if q == 2 # Two qubit MUBS 
        q1 = qubits[qubits_id[1]]
        q2 = qubits[qubits_id[2]]
        if no == 1
            circuit.h(q1)
            circuit.h(q2)
            circuit.s(q2)
            circuit.cx(q2,q1)
            circuit.s(q1)
            circuit.cx(q2,q1)
            return
        elseif no == 2
            circuit.h(q1)
            circuit.h(q2)
            circuit.s(q1)
            circuit.cx(q2,q1)
            circuit.s(q1)
            circuit.cx(q2,q1)
            return 
        elseif no ==3
           circuit.h(q1)
           circuit.h(q2)
           return
        elseif no == 4
            circuit.h(q1)
            circuit.h(q2)
            circuit.s(q1)
            circuit.s(q2)
            return
        else
            return # no circuit needed for Z
        end
    elseif q == 1 # One Qubits MUBS
        q1 = qubits[qubits_id[1]]
        if no == 1
            circuit.h(q1)
            return
        elseif no ==2
            circuit.h(q1)
            circuit.s(q1)
            return
        else
            return # Nothing needed for z
        end
    else
        @warn("passed in an invalid experiment, qubit number $q")
    end
            
end


"""
addReverseCircuit(circuit,experiment_type,qubits_id)

## Arguments
   circuit - the qiskit circuit we will add gates to.
   experiment-type: (Integer,Integer) - the stabilizer group we wont to generate (see below)
   qubits_id: Array of Int - the 'number' of the qubits the circuit is to be on.

Effectively this 'reverses' the circuit added by addCircuit. 

e.g. 

```
addCircuit(my_circuit,(2,1),[1,2])
addReverseCircuit(my_circuit,(2,1),[1,2,])
```

Will create a 'nothing' circuit i.e. self inverting.

The stabilizer circuits vary depending on the experiment-type. This is a tuple, the first number indicating the number of 
qubits (currently 1 or 2), the second which experiment currently 1:3 for single qubits and 1:5 for two qubits.

The qubits_id will be used to work out which qubits in the circuit to use. They are "JULIA INDEXED".

The qubits are identified via the circuit ```qubits = circuit.qubits```. But because of the way PyCall works, they will be 
in a 1 indexed array. e.g. qiskitits qubit_0 will be qubits[1].

Therefore if you want the circuit to be on qubits 0 and 1, you pass in [1,2] into qubits_id.

The circuits that will be added will generate the following for two qubits its:

1.    [[0,0,0,0],[1,1,0,1],[1,0,1,1],[0,1,1,0]], #II ZX YZ XY
2.    [[0,0,0,0],[1,1,1,0],[0,1,1,1],[1,0,0,1]], #II ZY XZ YX
3.    [[0,0,0,0],[0,0,0,1],[0,1,0,0],[0,1,0,1]], #II IX XI XX
4.    [[0,0,0,0],[0,0,1,0],[1,0,0,0],[1,0,1,0]], #II IY YI YY
5.    [[0,0,0,0],[0,0,1,1],[1,1,0,0],[1,1,1,1]]] #II IZ ZI ZZ

For one qubit its

1.    [[0,0],[0,1]], # IX
2.    [[0,0],[1,0]], # IY
3.    [[0,0],[1,1]], # IZ

## Returns nothing

The qiskit circuit will have had the gates added.
"""                 
function addReverseCircuit(circuit,experiment_type,qubits_id)
# Should probably actually reverse the circuit - this is 
# possibly simpler to understand but will be harder to maintain.
    qubits = circuit.qubits
    (q,no) = experiment_type
    if q == 2 # Two qubit MUBS 
#              [[0,0,0,0],[1,1,0,1],[1,0,1,1],[0,1,1,0]], #II ZX YZ XY
#              [[0,0,0,0],[1,1,1,0],[0,1,1,1],[1,0,0,1]], #II ZY XZ YX
#              [[0,0,0,0],[0,0,0,1],[0,1,0,0],[0,1,0,1]], #II IX XI XX
#              [[0,0,0,0],[0,0,1,0],[1,0,0,0],[1,0,1,0]], #II IY YI YY
#              [[0,0,0,0],[0,0,1,1],[1,1,0,0],[1,1,1,1]]] #II IZ ZI ZZ
        
        q1 = qubits[qubits_id[1]]
        q2 = qubits[qubits_id[2]]
        if no == 1
            circuit.cx(q2,q1)
            circuit.s(q1)
            circuit.s(q1)
            circuit.s(q1)
            circuit.cx(q2,q1)
            circuit.s(q2)
            circuit.s(q2)
            circuit.s(q2)
            circuit.h(q1)
            circuit.h(q2)
            return
        elseif no == 2
            circuit.cx(q2,q1)
            circuit.s(q1)
            circuit.s(q1)
            circuit.s(q1)
            circuit.cx(q2,q1)
            circuit.s(q1)
            circuit.s(q1)
            circuit.s(q1)
            circuit.h(q1)
            circuit.h(q2)
            return 
        elseif no ==3
           circuit.h(q1)
           circuit.h(q2)
           return
        elseif no == 4
            circuit.s(q1)
            circuit.s(q2)
            circuit.s(q1)
            circuit.s(q2)
            circuit.s(q1)
            circuit.s(q2)
            circuit.h(q1)
            circuit.h(q2)
            return
        else
            return # no circuit needed for Z
        end
    elseif q == 1 # One Qubits MUBS
        q1 = qubits[qubits_id[1]]
#                     [[0,0],[0,1]], # IX
#                     [[0,0],[1,0]], # IY
#                     [[0,0],[1,1]], # IZ
        if no == 1
            circuit.h(q1)
            return
        elseif no ==2
            circuit.s(q1)
            circuit.s(q1)
            circuit.s(q1)
            circuit.h(q1)
            return
        else
            return # Nothing needed for z
        end
    else
        @warn("passed in an invalid experiment, qubit number $q")
    end
            
end
   
"""
setUpExperiment(experiment,circuit;reverse=false)

## Arguments
    experiment: List of tuples - the stabilizers to be extracted
    circuit: circuit to be done
    reverse: apply or removing? Boolean.

Cycles through the experiment applying the correct gates in the circuit.
The experiment is set up as a list of tuples (q_no,exp_no), representing
the number of qubits (currently 1 or 2) and the type of experiment.

For example [(1,3),(2,1),(2,3)] would be a 5 qubit experiment.
-  The 1 qubit experiment 3 would be applied to the first qubit.
-  The 2 qubit experiment 1 would be applied to the second and third qubits
-  The 2 qubit experiment 3 would be applied to the fourth and fifth qubits.

## Returns

Nothing, but amends the supplied circuit.
"""
function setUpExperiment(experiment,circuit;reverse=false)
    numberOfQubits = foldl((a,b)->a[1]+b[1],experiment)
    start = 1
    for (q,no) in experiment
        if q ==2
            ids = [start,start+1]
            start +=2
            if reverse
                addReverseCircuit(circuit,(q,no),ids)
            else
                addCircuit(circuit,(q,no),ids)
            end
        elseif q==1
            ids = [start]
            start +=1
            if reverse
                addReverseCircuit(circuit,(q,no),ids)
            else
                addCircuit(circuit,(q,no),ids)
            end
        else
            @warn("Currently we only deal with local mubs, not on $q qubits")
        end
    end
end

"""
genPTwirl(circuit,len,qubit_id)

generates a Pauli twirl on the circuit consisting of len gates + inversion gate. On qubit_id

##Arguments
	circuit - the qiskit circuit to alter
	len - an integer representing the number of gates in the twirl
	qubit_id - an integer representing the qubit. Index from 1 (Julia), the qubit register is retrieved from circuit and is not needed.
## Returns
	Nothing - the circuit is altered.
"""
function genPTwirl(circuit,len,qubit_id)
    qubits = circuit.qubits
    list = rand(1:4,len)
    total = foldr(*,map(x->superPaulis[x],list))
    push!(list,findfirst(x->x==total,superPaulis))
    q1 = qubits[qubit_id]
    for i in list
        if i == 2
            circuit.x(q1)
        elseif i == 3
            circuit.y(q1)
        elseif i == 4
            circuit.z(q1)
        else
            circuit.id(q1)
        end
        circuit.barrier(q1)
    end

end

"""
getCircuit(experiment,noOfGates)

##Arguments
	experiment: a list of tuples (such as that obtained from generateSensibleSubsamples)
	noOfGates: the number of Paulis we want in the middle (to twirl)

Generates a qiskit circuit that implements the experiment passed in (e.g [(2,1) (2,2)] will generate
a cricuit that uses two 2 qubit MUBS at the beginning and end (1 and 2 respectively) and then performs
a Pauli Twirl in the middle). The circuits are self-inverting.

Appropriate barriers are inserted.

##Example
```
(experiments,paulisAll,ds) = generateSensibleSubsamples([(2,1),(1,2)]) # Get some experiments
circuit = getCircuit(experiments[1],12) # twelve pauli twirl circuit.
qiskit.execute(circuit,qiskit.Aer.get_backend("qasm_simulator"),shots=100).result().get_counts() # noisless simulator.
circuit.draw() # Will draw the circuit
```

## Returns
 	the qiskit circuit.
 """
function getCircuit(experiment,noOfGates)
	numberOfQubits = -1
	if length(experiment) == 1
		numberOfQubits = experiment[1][1]
	else
    	numberOfQubits = foldl((a,b)->a[1]+b[1],experiment)
    end
    circuit =qiskit.QuantumCircuit(numberOfQubits,numberOfQubits)
    setUpExperiment(experiment,circuit)
    circuit.barrier()
    for i in 1:numberOfQubits
        genPTwirl(circuit,noOfGates,i)
    end
    setUpExperiment(experiment,circuit,reverse=true)
    circuit.barrier()
    circuit.measure(0:numberOfQubits-1,0:numberOfQubits-1)
    return circuit
end

"""
doASeries(experiment,lengths,sequences,shots,noise_mode)
## Arguments
    experiment - the experiment to perform (see getCircuit for detailed explanation)
    lengths - a list of the lengths to perform the experiment on.
	sequences - the number of randomised sequences to perform at each length.
	shots - the number of measurements of each sequences
	noise-model - the noise model you want passed to the 'Aer' simulator.

A convenience function to simulate the basic information gathering algorithm.
We assume a an Aer simulator, and you supply the noise model.
For edifferent simulators or to run on a real device you will have to re-implement this function, probably
just change the execute line.
## Returns
	a list of the a list of the results of each sequence.(ordered by the same way as lengths).
Note on a real device you probably want to randomise the order that different lengths are called i.e. don't do all the short ones, slightly longer and then longer ones in a row.
"""
function doASeries(experiment,lengths,sequences,shots,noise_model)
    results = []
    for l in lengths
        perSeqResults = []
        for seq = 1:sequences
            circuit = getCircuit(experiment,l)
            counts = qiskit.execute(circuit,qiskit.Aer.get_backend("qasm_simulator"),noise_model=noise_model,shots=shots).result().get_counts()
            push!(perSeqResults,counts)
        end
        push!(results,perSeqResults)
    end
    return results
end

"""
getCounts(results,lengths,noOfQubtis)

Convenience funciton that takes the results from doASeries and consolidates the counts accross different sequences
i.e. if we have 10 sequences at lengths 10,20 and 40, groups all the sequences at length 10 together, the seqeunces at length 20 ... etc.
It then transforms into a probability vector (i.e. the vector for each length by the total number of counts).

NOTE: we assume that the machine is able to handle a vector of size 2^noOfQubits

##  Returns
	A list of probability vectors, one for each length.
"""
function getCounts(results,lengths,noOfQubits)
    consolidatedCounts = [zeros(2^noOfQubits) for _ in 1:length(lengths)];
    for (ix,seqs) in enumerate(results)
        for s in seqs
            for k in keys(s)
                consolidatedCounts[ix][parse(Int,k,base=2)+1] += s[k]
            end
        end
    end
    #Normalise to a percentage.
    for i  in 1:length(consolidatedCounts)
        consolidatedCounts[i]  = consolidatedCounts[i]./sum(consolidatedCounts[i])
    end
    return consolidatedCounts
end

"""
extractEigenvalues(counts,lengths)

Uses the functionality of Juqst - quantumNoise algorithms to extract the eigenvalues from a vector of counts.
For documentation on fitTheFidelities, see rharper2/Juqst.jl on gitHub.
To see how the counts should look see getCounts or the example workbooks.

## Returns 
	the Eigenvalues for the supplied probability vector.
"""
function extractEigenvalues(counts,lengths)
    (params,l, failed) = fitTheFidelities(lengths,counts)
    return vcat(1,[p[2] for p in params]) # We don't fit the first one, it is always 1 for CPTP maps
end



function binaryArrayToNumber(x)
    return foldl((y1,y2)->y1<<1+y2,x)
end

function generateEigenvalueRecur(stabs,current)
    if length(stabs) != 1
        current = generateEigenvalueRecur(stabs[2:end],current)
        q = stabs[1]
        multiplier = 2^length(q)
        newCurrent = []
        for c in current
            for i in [binaryArrayToNumber(x) for x in stabs[1]]
                push!(newCurrent,multiplier*c+i)
            end
        end
        return newCurrent
    else
        return [binaryArrayToNumber(x) for x in stabs[1]]
    end
end

"""
generateEigenvalues(stabs)

## Arguments
	stabs - a list of the Paulis used in the experiment, this can be obtained from the generateSensibleSubsamples function.

Uses a recursive algorithm to work out which global eigenvalues will be learnt through the supplied Paulis/experiment. You can use
PEEL.fidelityLabels(ix-1) to see the Pauli representation of these eigenvalues. Note the need to subtract 1, to mvoe from Julia 1 index
to a 0 index label (i.e. in Julia II...I is 1, whereas fidelityLabels expects this to be 0).
## Returns
	a list of eigenvalues (Julia 1 indexed).
"""
function generateEigenvalues(stabs)
    x = generateEigenvalueRecur(stabs,[])
    return x.+1 # Julia indexing.
end

"""
getStabiliserForExperiment(experiment)
	
For a given experiment (see generateSensibleSubsamples)  returns the Paulis (in bit form) that are accessed via the experiment.
"""
function getStabilizerForExperiment(experiment)
    paulis =[]
    for (q,s) in experiment
        if q == 2
            push!(paulis,all2QlMuBs[s])
        else
            push!(paulis,potentialSingles[s])
        end
    end
    return paulis
end



