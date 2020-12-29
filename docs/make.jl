using Documenter, FFTInterpolations 

 # set seed fixed for documentation

makedocs(modules = [FFTInterpolations], 
         sitename = "FFTInterpolations.jl", 
         pages = Any[
            "FFTInterpolations.jl" => "index.md",
         ]
        )

deploydocs(repo = "github.com/roflmaostc/FFTInterpolations.jl.git")
