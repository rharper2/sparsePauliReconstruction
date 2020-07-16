using Documenter

include("../peel.jl")
using Main.PEEL
include("../localPeelFunctions.jl")

makedocs(sitename="Scalable Estimation Documentation",
	    pages = [
        	"Introduction" => "index.md",
        	"Code" => Any[
            	"Main functions" => "peel.md",
            	"Local stabiliser functions" => "localPeel.md"
        	],
        	"Example Workbooks" => "example.md"
    	]
 )

