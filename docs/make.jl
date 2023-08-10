using QuadratureRules
using Documenter

makedocs(;
    modules=[QuadratureRules],
    authors="Michael Kraus",
    repo="https://github.com/juliagni/QuadratureRules.jl/blob/{commit}{path}#L{line}",
    sitename="QuadratureRules.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://juliagni.github.io/QuadratureRules.jl",
        assets=String[],
    ),
    pages=[
        "Home"    => "index.md",
        "Library" => "library.md",
    ],
)

deploydocs(;
    repo   = "github.com/JuliaGNI/QuadratureRules.jl",
    devurl = "latest",
    devbranch = "main",
)
