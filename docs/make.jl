using QuadratureRules
using Documenter
using DocumenterCitations

DocMeta.setdocmeta!(QuadratureRules, :DocTestSetup, :(using QuadratureRules); recursive=true)

bib = CitationBibliography(joinpath(@__DIR__, "src", "references.bib"))

makedocs(;
    plugins=[bib],
    modules=[QuadratureRules],
    authors="Michael Kraus",
    repo=Remotes.GitHub("JuliaGNI", "QuadratureRules.jl"),
    sitename="QuadratureRules.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://juliagni.github.io/QuadratureRules.jl",
        assets=String[],
    ),
    pages=[
        "Home"                 => "index.md",
        "Numerical Quadrature" => "quadrature.md",
        "Quadrature Rules"     => "rules.md",
        "Library"              => "library.md",
        "References"           => "references.md",
    ],
)

deploydocs(;
    repo   = "github.com/JuliaGNI/QuadratureRules.jl",
    devurl = "latest",
    devbranch = "main",
)
