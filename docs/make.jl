using FFTInterpolations, Documenter 

 # set seed fixed for documentation
DocMeta.setdocmeta!(FFTInterpolations, :DocTestSetup, :(using FFTInterpolations); recursive=true)
makedocs(modules = [FFTInterpolations], 
         sitename = "FFTInterpolations.jl", 
         pages = Any[
            "FFTInterpolations.jl" => "index.md",
         ]
        )

deploydocs(repo = "github.com/roflmaostc/FFTInterpolations.jl.git", devbranch="main")
