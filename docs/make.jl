using FFTResampling, Documenter 

 # set seed fixed for documentation
DocMeta.setdocmeta!(FFTResampling, :DocTestSetup, :(using FFTResampling); recursive=true)
makedocs(modules = [FFTResampling], 
         sitename = "FFTResampling.jl", 
         pages = Any[
            "FFTResampling.jl" => "index.md",
         ]
        )

deploydocs(repo = "github.com/roflmaostc/FFTResampling.jl.git", devbranch="main")
