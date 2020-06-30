# Copyright Robin Harper 2019-2020


module PEEL

using LsqFit, Hadamard, DelimitedFiles, LinearAlgebra
using PyPlot, Statistics, SparseArrays
import Base.show

# Defines a helper function in jupyter notebooks its /oplus<tab>
⊕ = (x,y)->mod.(x+y,2)

export generateFromPVecSamples4N,twoPattern,fourPattern,simplexUp,generatePauliFromClifford
export getχ,bitArrayToNumber





function getχ(a,b)
    return count_ones(a&b) %2 == 0 ? +1 : -1
end

function bitArrayToNumber(toS)
    return foldl((x,y)->x<<1+y,toS)
end



# The way I ended up doing it in the FULL MUB is probably better, but this is easy to follow

"""
fourPattern(p)

## Arguments
     ```p Array{Array{Array{Int64,1},1},1}``` represents a series of commuting two qubit paulis.

Given a MUB works out what binary string will go in what bin
This takes a 4 bit (two Pauli) set of commuting Paulis, eg.

```
1-element Array{Array{Array{Int64,1},1},1}:
 [[0, 0, 0, 0], [1, 0, 1, 1], [0, 1, 1, 0], [1, 1, 0, 1]]
 ```

And then returns the Array that shows what bin each of the 16 two qubit Paulis go into e.g.

```
[[["0000", "0111", "1001", "1110"], ["0001", "0110", "1000", "1111"], ["0011", "0100", "1010", "1101"], ["0010", "0101", "1011", "1100"]]]
```

Which tells us e.g. that "1001" goes into bin 0 and "1101" goes into bin 2.

Returns an Array, containing a single array of 4 vectors each of size 4 --- splitting all 16 Paulis (rep as binary strings) into 4 bins.
If you want to think of this as a hashing function, it shows how to take 4 bits -> 2 bits using the supplied pattern.
"""
function fourPattern(p)
    allPatterns=[]
    for x in p
        patterns =[Hadamard.hadamard(16)[x+1,:] for x in [foldl((x,y)->x<<1+y,t) for t in x]]
        push!(allPatterns,
            map(y->reverse.(string.(findall(x->x!=0,y).-1,base=2,pad=4)),Hadamard.hadamard(4)*patterns))
    end
    return allPatterns
end

"""
twoPattern(p)

### Arguments
     ```p Array{Array{Int64,1},1}``` represents a series of commuting two qubit paulis.

Given a MUB works out what binary string will go in what bin
This takes a 2 bit (one Pauli) set of commuting Paulis, eg.

```
1-element Array{Array{Array{Int64,1},1},1}:
 [0, 0], [1, 1]
 ```

And then returns a two element Array that shows what bin each of the 16 two qubit Paulis go into e.g.

```
["00", "11"]
["01", "10"]
```

Which tells us e.g. that "11" goes into bin 0 and "01" goes into bin 2.

Returns an Array, containing an array of 2 vectors each of size 2 --- splitting all 4 Paulis (rep as binary strings) into 2 bins.
If you want to think of this as a hashing function, it shows how to take 2 bits -> 1 bits using the supplied pattern.
"""
function twoPattern(ps)
    patterns =[Hadamard.hadamard(4)[x+1,:] for x in [foldl((x,y)->x<<1+y,t) for t in ps]]
    return map(y->reverse.(string.(findall(x->x!=0,y).-1,base=2,pad=2)),Hadamard.hadamard(2)*patterns)
end


""" 
getIndexOf(hadamardMap,no,dictMap=Dict())
### Arguments
	```hadamardMap: An array``` of hashing bins (see later)
	```no: Integer ``` the index of the Pauli we want to know which bin it gotes into 
	```dictMap: Dictionary``` optional mapping eg qubit 2->qubit 4.

 Given a map (constructed from fourPattern and twoPattern), this will
 split the qubits up into the pairs/singles identified by the map and
 then determine the hash destingation of a particular Pauli, i.e. which bin it goes into.
 The optional dictMapping, allows us to map qubits to different qubits. This might be 
 relevant for different topologies. E.g. we might want the pair to span qubits 1 and 14, rather
 than 1 and 2, so we would map 2 onto 14. There should be examples later.

 **NOTE** because we use the least significant (rh) qubit - there is a subte reversal of bits. See example.

### Example

The hadamard Map for a six qubit system might look like this
```
4-element Array{Any,1}:
 [["00", "10"], ["01", "11"]]
 [["0000", "0110", "1101", "1011"], ["1100", "1010", "0001", "0111"], ["1000", "1110", "0101", "0011"], ["0100", "0010", "1001", "1111"]]
 [["0000", "1110", "1001", "0111"], ["1100", "0010", "0101", "1011"], ["0100", "1010", "1101", "0011"], ["1000", "0110", "0001", "1111"]]
 [["00", "01"], ["10", "11"]]
 ```

 Where we have used twoPattern for qubits 0 and 5 and fourPattern for qubits 1&2 and 3&4

 getIndexOf(patternAbove, 23), would then return 6.

 Because 23 = "000000010111"
 Then we map the bins as 
 00-> reversed to 00 -> bin 0 -> 0
 0000 -> reversed to 0000 ->  bin 0 -> 00
 0101 -> reversed to 1010 ->  bin 2 -> 10
 11 ->  reversed to 11 -> bin 1 -> 01

 Therefore the index is 0 00 10 01 = 5, as we one index we return 6.

 ## Returns

 Index of bin (1 indexed).

"""
function getIndexOf(hadamardMap,no,dictMapping=Dict())
    bitsize = sum([length(x)>1 ? length(x) : floor(Int,log2(x[1])) for x in hadamardMap]) # 2* number of qubits
    bitst = (string(no,base=2,pad=bitsize))
    if length(dictMapping)!=0
        currentValues = [parse(Int,x) for x in bitst]
        newB = copy(currentValues)
        for i = 1:length(dictMapping)
               currentValues[dictMapping[i]*2-1] = newB[i*2-1]
               currentValues[dictMapping[i]*2] = newB[i*2]
        end
        bitst = string(foldl((x,y)->x<<1+y,currentValues),base=2,pad=bitsize)
    end
    bitsize = sum([length(x)>1 ? length(x) : floor(Int,log2(x[1])) for x in hadamardMap]) # 2* number of qubits
    #lshift = round(Int,bitsize/2)-1
    occupancy = [length(x)>1 ? length(x) : x[1] for x in hadamardMap]
    cpos = 1
    indexArray = []
    for hm in hadamardMap
        if (length(hm)>1) # We assume that we have a list of MUBS
            bits = bitst[cpos:cpos+length(hm)-1]
            mv = findfirst(x->x==true,in.(reverse(bits),hm))-1
            push!(indexArray,mv)
            cpos = cpos+length(hm)
        else # we have a binary watchamacallit
            numberOfBits = floor(Int,log2(hm[1]))
            bits = bitst[cpos:cpos+numberOfBits-1]
            push!(indexArray,parse(Int,bits,base=2))
            cpos = cpos+numberOfBits
        end
                
    end

    #print(indexArray)

    val = indexArray[1]
    hmLengths = [floor(Int,log2(t)) for t in occupancy]
    st = ""
    for ix in 1:length(hmLengths)
        st = st*string(indexArray[ix],base=2,pad=hmLengths[ix])
    end
    return parse(Int,st,base=2)+1
end


"""
constructTheKeyIndexMUB(indx,listOfPs,map=[])

## Arguments
    ```indx``` The index of a Pauli, we want to ascertain which bins that Pauli was hashed into.
	```listOfPs``` The used to sample/bin the Paulis
	```map``` A map of qubit substitutions e.g mapping qubit 2 to 14

Given a Paulis to test, returns the index into all the bins that Paulis is hashed into 
(it is assumed we have a number a sub-sampling groups).
Basically just iterates over them. If we supply a map, it should be an
array of maps, one for each different set of n-qubit MUBS supplied.

It uses getIndexOf to do all the heavy lifting, the documentation for that function 
provides an example of the type of map it requires, this is an array of these maps.

Returns and array of offsets representing the index the Pauli gets hashed to.

"""
function constructTheKeyIndexMUB(indx,listOfPs;map=[])
    if map == []
        return [getIndexOf(x,indx) for x in listOfPs]
    else
        return [getIndexOf(x,indx,map[ix]) for (ix,x) in enumerate(listOfPs)]
    end
end


""" 
generateFromPVecSamples4N(pvec,d=[],dictMapping)

## Arguments
	```pvec``` The list of Paulis which are used in our experiment. Given as an array of stabilisers (see example)
	```d``` The offset, this is the binary string 'hashed' into the Paulis given by pvec.
	```map``` an optional dictionary mapping qubits to other qubits.

Generates the fidelities we are going to sample from for a given offset 
vector (d).

By way of example, if we supplied as a pvec

```
[[0, 0, 0, 0], [0, 1, 1, 1], [1, 0, 0, 1], [1, 1, 1, 0]]

```
 Then we would get 0,7,9 and 14.

 If we also suplied a d, [1 0 0 0], then we get 8, 15, 1 6

 A pvec as 

 ```
 [[0, 0, 0, 0], [0, 1, 1, 1], [1, 0, 0, 1], [1, 1, 1, 0]]
 [[0, 0, 0, 0], [0, 1, 1, 1], [1, 0, 0, 1], [1, 1, 1, 0]]
 ```

 Would give us 0, 7, 9, 14, 112, 119, 121, 126, 144, 151, 153, 158, 224, 231, 233, 238

 ### Returns
 	Array{Int64,1} with the relevant fidelity indexes. (zero indexed).

"""
function generateFromPVecSamples4N(pvec,d=[],dictMapping=Dict())
    subSamples = []
    if d == []
        d = zeros(Int64,sum([length(x[1]) for x in pvec]))
    end
    
    maxValue = foldl(*,[length(x) for x in pvec])
    occupancy = vcat(map(x->length(x),pvec),1)
    # Tuple of what we need to divide by and what we need to mod by
    # So for example if the lengths given were [4,8,4]
    # Then we need to to range from 0:(4*8*4=256)-1
    # We need to divide our 'c' by (8*4)%4 and (4)%8 and (1)%4 to get the three indexes needed.
    divs = [(foldl(*,occupancy[ix:end]),occupancy[ix-1]) for ix=length(occupancy):-1:2]
    for c = 0:maxValue-1
        # express c as a four  array of length 1
        iv = reverse(map(x->floor(Int,c/x[1])%x[2]+1,divs))
#           print("$iv\n")
        xs_o = foldl(vcat,map(x->pvec[x[1]][x[2]],enumerate(iv)))
        xs = copy(xs_o)
        if length(dictMapping) != 0
           for i = 1:length(dictMapping)
               xs[dictMapping[i]*2-1] = xs_o[i*2-1]
               xs[dictMapping[i]*2] = xs_o[i*2]
           end
        else
                xs = xs_o
        end
            #println(xs)
        nxs = xs ⊕ d
            #print("$(itoP(foldl((x,y)->x<<1+y,nxs)))\n")
        push!(subSamples,foldl((x,y)->x<<1+y,nxs))
    end
    return [x for x in subSamples]
end


"""
closeToZero(bin,BITS,cutoff=1e-15)

## Arguments
	```bin``` An array of values in the bin and its various offsets. 
	```bits``` The number of bits (typically 2*#ofqubits)
	```cutoff``` the entropy in the bin, below which we will assume it is zero.

Checks if the values in the bin have an absolute value close to zero, where we sum the 
squares of the entries in the bin, divide by (BITS+1) and compare to the cutoff.

## Returns
	true or false.
"""
function closeToZero(indexes,BITS;cutoff=1e-15)
    # checks to see if the supplied [U[k̂] for randomDs] are close to zero. i.e. no weights
    # Just now , this is effectively ||U||^2 
    U = [ix^2 for ix in indexes]
    if (sum(U)/(BITS+1) < cutoff)
        return true
    else
        return false
    end
end

""" 
checkAndExtractSingleton(indexes,BITS;cutoff=0.000002)

## Arguments
	```indexes``` A list of arrays, being the bins and their offsets.
	```Bits``` the number of bits (being 2*qubits)
	```cutoff``` what the entropy has to be less than
Checks the 'entropy' in the bin (ie sum of the square of the fluctations around the mean of the abs value)/BITS
If it is less than the cutoff admits it as a singleton.
Extracts the Pauli (using majority vote if given multiple samples)
-  Note if you wanted to implement the "Error Correcting Code in the offsets" as described in the appendix of the paper then this funciton would be different.

## Returns
    A tuple of: (found,index, value), where found is boolean, the index is an integer and value is float (between 0 and 1).
  
"""
function checkAndExtractSingleton(indexes,BITS;cutoff=0.000002)
    stringsFound = zeros(BITS)
    threshold = floor((length(indexes)-1)/2)
    for (ii,i) in enumerate(indexes)
        #print("Check and Extract $ii : $i\n")
        pcheck = i[1]
        toAdd = ""
        for (ix,xx) in enumerate(i[2:end])
            if pcheck >=0
                if xx >=0 
                    stringsFound[ix]+=1
                end
            else
                if xx < 0 
                    stringsFound[ix]+=1
                end
            end
        end
    end
    #print("$stringsFound\n")
    kstr =  [x > threshold ? "0" : "1" for x in stringsFound]
    kii = foldl(*,kstr)
    kval = parse(Int,kii,base=2)
    # so our estimate of this particular string should be indexes[1][1]
    Uk̂ = indexes[1][1]
    Uk̂bar = mean(mean([abs.(x) for x in indexes]))
    entropy = sum([sum([(abs.(x)-Uk̂)^2 for x in ix])/(BITS+1) for ix in indexes])
    #print(entropy)
    if entropy < cutoff
       # print("Singleton $(Uk̂) and $(mean(mean([abs.(x) for x in extracted])))")
        return (true,kii,Uk̂bar)
    else
        return (false,0,0)
    end

end

"""
peelBack(listOfX,listOfPs,singletonBits,singletonValue,found,ds,mappings)

## Arguments
	```listOfX``` - the bins and their offsets <- this gets modified
	```listOfPs``` - the stabiliser groups used to create the bins (see generateFromPVecSamples4N)
	```singletonBits``` - the bit index of the Pauli to peel back
	```singletonVlaue``` - the value of the Pauli to peel back
	```found``` - a Dict of Pauli bitstrings that have been found and their value
	```ds``` - the offsets used.
	```mappings``` - an optional mapping of qubit->qubit.

Peel back algorithm tailored for noise version.

Check if we have already found the supplied Pauli (is it in found?)
    - If so check if the value we think we have found is greater than the previous sample
       - If not, then just return - we had already found it, and the noise has just confused us!
       - If yes, then the one we found was probably wrong and a result of the noise, just update the value.

*[Query this logic, it seems to work but perhaps, we may want to update the bins with the difference.]*

    - Otherwise (i.e. we hadn't found it before.)
    	- Get the index of that Pauli into each of the other bin sets
		- And add or subtract (depending on relevant offset) the supplied value of the Pauli
		- Update found by adding in this Pauli and the value.
"""
function peelBack(listOfX,listOfPs,singletonBits,singletonValue,found,ds,mappings)
    if (singletonBits in keys(found))
        #print("We already done that one!\n")
        vval = parse(Int,singletonBits,base=2)
        vstring = singletonBits
        if singletonValue > mean(found[vstring])
                found[vstring] = [singletonValue]
        end
        return
        #@assert false
    end
    vval = parse(Int,singletonBits,base=2)
    vstring = singletonBits
    found[vstring] = vcat(get(found,vstring,[]),[singletonValue])
    # Get the keys/bin number for all the other hash/MUBs we have.
    theKeys = constructTheKeyIndexMUB(vval,listOfPs,map=mappings)
    for dictIndex in 1:length(theKeys)
        pv = theKeys[dictIndex]
        #print("\t****** Stepping through left to right nodes (so to speak) for identified singleton\n")
        extractValue = pv#parse(Int,pv,base=2)+1
        extracted = []
            push!(extracted,[x[extractValue] for x in listOfX[dictIndex]])
        newIndex = []
        for (sidx,bd) in enumerate(ds)
               push!(newIndex,extracted[1][sidx] + -1*PEEL.getχ(PEEL.bitArrayToNumber(bd),vval)*singletonValue)
        end
         theX = listOfX[dictIndex]
         for (ix1) in 1:length(theX)
             theX[ix1][extractValue] = newIndex[ix1]
         end
    end 
end

# ======================================================================================
# The following are just convience functions for giving us labels we can use
# ======================================================================================


"""
probabilityLabels(index; qubits=2)

## Arguments
    ```index (Int)``` The index (0 index) for which we want the label.
    ```qubits (Int)``` The number of qubits in the label, defaults to 2
Because we use the Walsh-Hadmard transform based on bits (i.e. the standard 'natural' transform rather than a Pauli-commuting transform)
The transform changes the way we need to label our qubits. In the probability regime we have the following bit labelling:
```
	00-> I
	01-> Y
	10-> X
	11-> Z
```
The first qubit is the least signifcant number (i.e. in binary would appear on the right hand side).
So the number 1001 is XY, where qubit 0 = Y.

### Returns 
The probability label for the index.

"""
function probabilityLabels(x;qubits=2)
    str = string(x,base=4,pad=qubits)
    paulis = ['I','Y','X','Z']
    return map(x->paulis[parse(Int,x)+1],str)
end

"""
fidelityLabels(index; qubits=2)

## Arguments
    ```index (Int)``` The index (0 index) for which we want the label.
    ```qubits (Int)``` The number of qubits in the label, defaults to 2
Because we use the Walsh-Hadmard transform based on bits (i.e. the standard 'natural' transform rather than a Pauli-commuting transform)
The transform changes the way we need to label our qubits. In the probability regime we have the following bit labelling:
```
	00-> I
	01-> X
	10-> Y
	11-> Z
```
The first qubit is the least signifcant number (i.e. in binary would appear on the right hand side).
So the number 1001 is YX, where qubit 0 = X.

### Returns 
The fideltiy label for the index.

"""
function fidelityLabels(x;qubits=2)
    str = string(x,base=4,pad=qubits)
    paulis = ['I','X','Y','Z']
    return map(x->paulis[parse(Int,x)+1],str)
end



# ======================================================================================
# The following functions are to generate a Pauli probability from a "Clifford" probability
# ======================================================================================


"""
simplexUp(n)

## Arguments

   `n (Int)`  = number to project onto 

Projects onto the n-simplex. i.e. returns n numbers such that the sum of them is 1.

## Returns

vector of n-numbers that add up to 1.

"""
function simplexUp(n)
    return diff(vcat(0,sort(rand(n-1)),1))
end


""" 
noOfOnesIn(ix)

## Arguments

	`ix (Int)` A number representing an index in Julia 

Return the number of '1' in the binary form of the number -1. 
Note that we subtract 1 from the supplied index to take into account Juila's 1 indexing. 
The 'index' 1 actually represents the number 0, which has 0 '1's in its binary form.
The indexes 2 and 3, both have 1 index in their binary form (e.g 01 and 10), index 4
(which is binary 11) has 2 and so on.

	Return the number of 1s in the binary representation of ix-1.
"""
function noOfOnesIn(ix)
    return length(filter(x->x=='1',string(ix-1,base=2)))      
end


"""
gen(n)

## Arguments 
     `n (Int) ` Number of qubits

Calculates and returns all possible (non-identity) Paulis  for a certain system size. The Paulis
are represented by the numbers 1,2 and 3.


### Returns 

A vector of the Paulis. 
"""
function gen(n)
    if n == 1 
        return [1,2,3]
    else 
        x = gen(n-1)
        return vcat(map(y->vcat([1],y),x),map(y->vcat([2],y),x),map(y->vcat([3],y),x))
    end
end


"""
splitUp(n)

## Arguments
    `n Integer`: An 'index' into the Cliffo
	
Splits the "Clifford average" into the Paulis that were averaged into it.
The starting point here is that if we assume we had single-Clifford twirls on an n-qubit machine
We have our noise supoeroperator, diagonalised into a superoperator matrix, Pauli basis, with 2^n distinct
entries (4^n entries in total).  If we imagine we labelled these distinct entries sequentially then the first entry would represent the 'all I' eigenvalue.
(These distinct entries are in fact the eigenvalues we recover from a Hadamard transform of the measurement outcomes, where we measure in the computational basis).
Next there would then be three averaged eigenvalues representing the I...IX, I...IY, I....IZ eigenvalues and so on. What we are interested in is
splitting that up to retrieve the eigenvalues that were averaged to give that number.

For example, the no-error starts IIIII.IIII = index 1
	
-	splitUp(1) = 1
-	splitUp(2) (Z) = 2,3,4 i.e. IX, IY, IZ
-	splitUp(3)  (ZI) = 5,9,13 ie. XI, YI, ZI
-	splitUp(4) (ZZ) = 6,7,8,10,11,12,14,15,16 i.e. XX,XY, XZ, YX,YY,YZ, ZX, ZY, ZZ

## Returns 
A vector of the Paulis that would have been averaged together to give that specific Clifford.
"""
function splitUp(n)
    if n == 1 
        return [1]
    end
    vals = []
    for subs in gen(noOfOnesIn(n))
        original = string(n-1,base=2)
        rec = zeros(Int64,length(original))
        count = 1
        for x in subs
            while original[count] != '1'
                rec[count]= 0
                count += 1
            end
            rec[count] = x
            count+=1
        end
        push!(vals,1+foldl((x,y)->x*4+y,rec))
    end
    return vals
end


"""
generatePauliFromClifford(n)

## Arguments
    `cliffordPs`: Array{Float64,1}

	
Splits the "Clifford average" which we recover from the protocol in Efficient Learning of quantum Noise arXiv:1907.13022. into 
a fake Pauli distribution, that would average to the same Clifford average. We are going to use this to test the algorithm.

We do this by realising that each member of our 2^n distribution was in fact an average of certain Paulis. We work out how many
Paulis went into the distribution (```splitUp(n)```) and then project the number onto an appropriate simplex (```simpleUp```).

This allows us to recreate a 4^n distribution that would average down to the supplied distribution. Of course there is a random 
element with the projection up stage, but it is *experimentally* inspired if the original averaged protocol came from an experiment.

## Returns 
A Pauli probability distribution that result in the supplied probability distribution if they were the errors caused by the 
avarged noise in the machine and that machine was subject to single qubit Clifford twirls.
"""        
function generatePauliFromClifford(cliffordPs)
    s = Dict{Int64,Float64}()
    nos=0
    totalToSub = 0
    for (ix,x) in enumerate(cliffordPs)
        if x>0
            split  = splitUp(ix)
            # Number of of paulis without a Z
            nsplit = 2^round(Int,(log(3,length(split))))
            increaseRatio = (length(split)+nsplit/2)/length(split)
            totalToSub += x*(increaseRatio-1)
            vals = simplexUp(length(splitUp(ix))).*(x*increaseRatio)
            for (jx,j) in enumerate(splitUp(ix))
                s[j]=vals[jx]
            end
        end
    end
    s[1] -= totalToSub
    lastValue = length(cliffordPs).^2
    s[lastValue] = get(s,lastValue,0)
    sp = sparsevec(s)
    return Vector(sp)
end


end # of module.







