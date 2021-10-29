pushfirst!(LOAD_PATH, joinpath(@__DIR__, ".."))

using Documenter, GCMMesh

makedocs(
    sitename="GCMMesh.jl",
    doctest=false,
    strict=false,
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        mathengine = MathJax(Dict(
            :TeX => Dict(
                :equationNumbers => Dict(:autoNumber => "AMS"),
                :Macros => Dict(),
            ),
        )),
    ),
    clean = true,
    modules=GCMMesh,
    pages = Any[
        "Home"=>"index.md",
        "Graphs"=>"Graphs.md",
        "Mesh"=>"Mesh.md",
        "Partition"=>"Partition.md",
    ],
)

deploydocs(repo = "github.com/CliMA/GCMMesh.jl.git", target = "build")
