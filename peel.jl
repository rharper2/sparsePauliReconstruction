module PEEL

using LsqFit, Hadamard, DelimitedFiles, LinearAlgebra
using PyPlot, Statistics, SparseArrays
import Base.show


⊕ = (x,y)->mod.(x+y,2)
export generateFromPVecSamples4N,twoPattern,fourPattern,itoP,simplexUp,generatePauliFromClifford
export getχ,bitArrayToNumber




# The following functions used to generate a Pauli probability from a "Clifford" probability


"""
Projects onto the n-simplex i.e. returns n
numbers such that the sum of them is 1.
"""
function simplexUp(n)
    return diff(vcat(0,sort(rand(n-1)),1))
end


""" 
	Return the number of 1s in the binary representation
"""
function noOfOnesIn(ix)
    return length(filter(x->x=='1',string(ix-1,base=2)))      
end


"""
Returns the Paulis (as 1,2,3) for a system of size n
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
	Splits the "cliffords" into Paulis start IIIII.IIII = index 1
	so splitUp(1) = 1
	splitUp(2) (Z) = 2,3,4 i.e. IX, IY, IZ
	splitUp(3)  (ZI) = 5,9,13 ie. XI, YI, ZI
	splitUp(4) (ZZ) = 6,7,8,10,11,12,14,15,16 i.e. XX,XY, XZ, YX,YY,YZ, ZX, ZY, ZZ
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

# ===========================================
# The following are just convience functions for converting between the various forms we use
# ===========================================

itop = Dict('0'=>'I','1'=>'Z','2'=>'X','3'=>'Y')
"""
 btoP(x)
 converts a binary string into the Pauli representation of it.
 Here (unlike the symplex/stabiliser representation) this is a linear binary string
 where each Pauli occupies two bits. We map as:

`    00->I
     10->X
     11->Y
     01->Z `


 the mapping symplifies any conversion we may want to do into symplex form.
"""
function btoP(x,qubits=14)
    return map(x->itop[x],string(parse(Int64,x,base=2),base=4,pad=qubits))
end

"""
    itoP(x)

    x: 
    qubit: specify the number of qubits to represent, default 14 (historical)

this is a convenience function similar to btoP - except that it takes an integer rather than
a binary string. It is assumed the integer is formed from a binary string, where each
Pauli is mapped as:

```
    00->I
    10->X
    11->Y
    01->Z 
```


 then the binary string is coverted to base 10 integer. 
"""
function itoP(x,qubits=14)
    return map(x->itop[x],reverse(string(x,base=4,pad=qubits)))
end



function getχ(a,b)
    return count_ones(a&b) %2 == 0 ? +1 : -1
end

function bitArrayToNumber(toS)
    return foldl((x,y)->x<<1+y,toS)
end



# The way I ended up doing it in the FULL MUB is probably better, but this is easy to follow

"""
Given a MUB works out what binary string will go in what bucket
This takes a 4 bit (two Pauli) MUB.
As a hashing funciton the MUB takes 4 bits -> 2 bits.
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
Given a MUB works out what binary string will go in what bucket
This takes a 2 bit (one Pauli) MUB.
As a hashing funciton the MUB takes 2 bits -> 1 bits.
"""
function twoPattern(ps)
    patterns =[Hadamard.hadamard(4)[x+1,:] for x in [foldl((x,y)->x<<1+y,t) for t in ps]]
    return map(y->reverse.(string.(findall(x->x!=0,y).-1,base=2,pad=2)),Hadamard.hadamard(2)*patterns)
end

""" 
 Given a map (constructed from fourPattern and Two pttern), this will
 split the qubits up into the pairs/singles identified bu the map and
 then determine the hash funciton of a particular Pauli, i.e. which bucket it goes into.
 The optional dictMapping, allows us to map qubits to different qubits. This might be 
 relevant for different topologies. E.g. we might want the pair to span qubits 1 and 14, rather
 than 1 and 2, so we would map 2 onto 14. There should be examples later.
"""
function getIndexOf(hadamardMap,no,dictMapping=Dict())
    bitsize = sum([length(x)>1 ? length(x) : floor(Int,log2(x[1])) for x in hadamardMap]) # 2* number of qubits
    bitstring = (string(no,base=2,pad=bitsize))
    if length(dictMapping)!=0
        currentValues = [parse(Int,x) for x in bitstring]
        newB = copy(currentValues)
        for i = 1:length(dictMapping)
               currentValues[dictMapping[i]*2-1] = newB[i*2-1]
               currentValues[dictMapping[i]*2] = newB[i*2]
        end
        bitstring = string(foldl((x,y)->x<<1+y,currentValues),base=2,pad=bitsize)
    end
    bitsize = sum([length(x)>1 ? length(x) : floor(Int,log2(x[1])) for x in hadamardMap]) # 2* number of qubits
    #lshift = round(Int,bitsize/2)-1
    occupancy = [length(x)>1 ? length(x) : x[1] for x in hadamardMap]
    cpos = 1
    indexArray = []
    for hm in hadamardMap
        if (length(hm)>1) # We assume that we have a list of MUBS
            bits = bitstring[cpos:cpos+length(hm)-1]
            mv = findfirst(x->x==true,in.(reverse(bits),hm))-1
            push!(indexArray,mv)
            cpos = cpos+length(hm)
        else # we have a binary watchamacallit
            numberOfBits = floor(Int,log2(hm[1]))
            bits = bitstring[cpos:cpos+numberOfBits-1]
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
    Given a set of paulis to test, returns the index into all the buckets
    Basically just iterates over them. If we supply a map, it should be an
    array of maps, one for each different set of n-qubit MUBS supplied.
"""
function constructTheKeyIndexMUB(indx,listOfPs;map=[])
    #print("CTKIM $indx\n")
    if map == []
        return [getIndexOf(x,indx) for x in listOfPs]
    else
        return [getIndexOf(x,indx,map[ix]) for (ix,x) in enumerate(listOfPs)]
    end
end


""" 
    Generates the fidelities we are going to sample from for a given offset 
    vector (d).
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






function subsample(u,M,d=[])
    subSampleSize = size(M)[2]
    subSamples=[]
    if d == []
        d = zeros(size(M)[1])
    end
    
    for i in 1:2^subSampleSize
        # for some reason the paper has left as the lsb
        newi = map(x->parse(Int,x),collect(reverse(string(i-1,base=2,pad=subSampleSize))))
        toS = M*newi ⊕ d
        #print("$(string(foldl((x,y)->x<<1+y,(toS)),base=2,pad=14))\n")
        push!(subSamples,foldl((x,y)->x<<1+y,toS))
    end
    return [u[x+1] for x in subSamples]
end



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
    Checks the 'entropy' in the bucket
    and if it is less than the cutoff admits it as a singleton.
    Extracts the Pauli (using majority vote if given multiple samples)
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
    Peel back algorithm tailored for noise version.
    Check if we have already found it. (is it in found)
    If so check if the value we think we have found is greater than the previous sample
    If not, then just return - we had already found it, and the noise has just confused us!
    If yes, then the one we found was probably wrong and a result of the noise, just update the value.
            Query, we may want to update the buckets with the difference.
    Otherwise - we hadn't found it before. Get the index into each of the other bucket sets
    And add or subtract (depending on relevant offset)
    Update found.
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
    # Get the keys/bucket number for all the other hash/MUBs we have.
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



end







