using QuadratureRules
using Documenter
using DocumenterCitations

DocMeta.setdocmeta!(QuadratureRules, :DocTestSetup, :(using QuadratureRules); recursive = true)

bib = CitationBibliography(joinpath(@__DIR__, "src", "references.bib"))

makedocs(;
    plugins = [bib],
    modules = [QuadratureRules],
    authors = "Michael Kraus",
    repo = Remotes.GitHub("JuliaGNI", "QuadratureRules.jl"),
    sitename = "QuadratureRules.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://juliagni.github.io/QuadratureRules.jl",
        assets = String[],
        # the library reference gathers every docstring in the package on one page and so
        # exceeds the 100 KiB the default warning threshold allows, at 104 KiB as of this
        # writing, well short of the 200 KiB limit
        size_threshold_warn = 150 * 1024
    ),
    pages = [
        "Home" => "index.md",
        "Numerical Quadrature" => "quadrature.md",
        "Quadrature Rules" => "rules.md",
        "Library" => "library.md",
        "References" => "references.md"
    ]
)

deploydocs(;
    repo = "github.com/JuliaGNI/QuadratureRules.jl",
    devurl = "latest",
    devbranch = "main"
)
