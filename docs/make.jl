using Documenter
using CellBase

makedocs(
    sitename = "CellBase",
    format = Documenter.HTML(),
    modules = [CellBase]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
