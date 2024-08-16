using Documenter, QuantumGraining

#ENV["GKSwstype"] = "100" # enable headless mode for GR to suppress warnings when plotting

DocMeta.setdocmeta!(QuantumGraining, :DocTestSetup, :(using QuantumGraining); recursive = true)

const page_rename = Dict("developer.md" => "Developer docs") # Without the numbers

function nice_name(file)
  file = replace(file, r"^[0-9]*-" => "")
  if haskey(page_rename, file)
    return page_rename[file]
  end
  return splitext(file)[1] |> x -> replace(x, "-" => " ") |> titlecase
end

makedocs(;
  modules = [QuantumGraining],
  doctest = true,
  linkcheck = false, # Rely on Lint.yml/lychee for the links
  authors = "Leon Bello <lionbello@gmail.com> and contributors",
  repo = "https://github.com/leonbello/QuantumGraining.jl",
  sitename = "QuantumGraining.jl",
  format = Documenter.HTML(;
    prettyurls = true,
    canonical = "https://leonbello.github.io/QuantumGraining.jl",
    assets = ["assets/style.css"],
  ),
  pages = [
    "Home" => "index.md"
    "API" => "api.md"
    [
      nice_name(file) => file for
      file in readdir(joinpath(@__DIR__, "src")) if file != "index.md" && file != "api.md" && splitext(file)[2] == ".md"
    ]
  ],
)

deploydocs(; repo = "github.com/leonbello/QuantumGraining.jl", push_preview = true)

# makedocs(sitename="QuantumGraining.jl")
# pages = [
#         "index.md",
#         "theory.md",
#         "tutorial.md",
#         "api.md",
#         "reference.md",
#         "Examples" => [
#             ]
#     ]

#     makedocs(;
#     sitename = "QuantumGraining.jl",
#     modules = [QuantumGraining],
#     pages = pages,
#     checkdocs = :exports,
#     format = Documenter.HTML(
#                             mathengine=MathJax(),
#                             footer="[**Back to GitHub**](https://github.com/leonbello/QuantumGraining.jl)"
#                             )
#     )

# deploydocs(
#     repo = "github.com/leonbello/QuantumGraining.jl",
#     push_preview = false,
#     )
