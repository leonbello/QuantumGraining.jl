using Documenter, QuantumGraining

makedocs(sitename="QuantumGraining.jl")

ENV["GKSwstype"] = "100" # enable headless mode for GR to suppress warnings when plotting

pages = [
        "index.md",
        "theory.md",
        "tutorial.md",
        "api.md",
        "Examples" => [
            ]
    ]

    makedocs(
    sitename = "QuantumGraining.jl",
    modules = [QuantumGraining],
    pages = pages,
    checkdocs=:exports,
    format = Documenter.HTML(
                            mathengine=MathJax(),
                            footer="[**Back to GitHub**](https://github.com/leonbello/QuantumGraining.jl)"
                            )
    )

deploydocs(
    repo = "github.com/leonbello/QuantumGraining.jl.git",
    push_preview = false,
    )
