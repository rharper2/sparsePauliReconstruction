# Additional local Peel Functions

You need to include "localPeelFunctions.jl" to access these.
Assumes qiskit is enabled in your Julia set up.

I believe the easiest way to make qiskit available is as follows:

### To get qiskit - under Julia

You will need to install Conda, PyCall using the julia package manager.

Easiest way is the repl (hit ], then add Conda, add PyCall)
or...

```
using Pkg
Pkg.add("Conda")
Pkg.add("PyCall")
```

then use Conda to install pip
```
using Conda
Conda.add("pip")
```

Once that is done, then we can use pip to install qiskit.


```
using PyCall
run(`$(PyCall.pyprogramname) -m pip install qiskit`)
```

Then we can use qiskit from Julia!


## Main user functions

```@docs
generateSensibleSubsamples
```

```@docs
getCircuit
```


```@docs
doASeries
```

```@docs
getCounts
```


```@docs
extractEigenvalues
```

```@docs
generateEigenvalues
```


## More down in the weeds

```@docs
addCircuit
```

```@docs
addReverseCircuit
```

```@docs
setUpExperiment
```

```@docs
genPTwirl
```

```@docs
getStabilizerForExperiment
```
