using HermiteNormalForm
using Documenter

DocMeta.setdocmeta!(
    HermiteNormalForm,
    :DocTestSetup,
    :(using HermiteNormalForm);
    recursive = true,
)

makedocs(;
    modules = [HermiteNormalForm],
    authors = "Yingbo Ma <mayingbo5@gmail.com> and contributors",
    repo = "https://github.com/YingboMa/HermiteNormalForm.jl/blob/{commit}{path}#{line}",
    sitename = "HermiteNormalForm.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://YingboMa.github.io/HermiteNormalForm.jl",
        assets = String[],
    ),
    pages = ["Home" => "index.md"],
)

deploydocs(; repo = "github.com/YingboMa/HermiteNormalForm.jl", devbranch = "master")
