using DQGMRES
using Documenter

DocMeta.setdocmeta!(DQGMRES, :DocTestSetup, :(using DQGMRES); recursive=true)

makedocs(;
    modules=[DQGMRES],
    authors="Yonatan Delelegn",
    repo="https://github.com/yonatanwesen/DQGMRES.jl/blob/{commit}{path}#{line}",
    sitename="DQGMRES.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://yonatanwesen.github.io/DQGMRES.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/yonatanwesen/DQGMRES.jl",
    devbranch="main",
)
