
module SheapIO

export SheapOptions, SheapMetadata, run_sheap

using Printf
using Parameters

"""
Commandline options for SHEAP
"""
@with_kw struct SheapOptions
    format_in::String = "vec"
    metadata_file = nothing

    verbose::Bool = false
    quiet::Bool = false
    dim::Int = 3
    scale::Bool = true
    similarity_threshold::Float64 = 0.0
    cost_threshold::Float64 = 0.0
    use_tsne::Bool = false
    use_umap::Bool = false
    perplexity::Float64 = 15.0
    knn::Int = 15
    uniform_random::Bool = false
    uniform_random_packing::Float64 = 0.5
    compression_factor::Float64 = 10.0
    kl_divergence_loss::Bool = false
    cross_entropy_loss::Bool = true
    hard_sphere::Bool = false
    hard_sphere_core_strength::Float64 = 1.0
    hard_sphere_steps::Int = 500
    hard_sphere_tol::Int = 2.0
    sphere_radius::Float64 = 0.02
end

"""
Record the default options of SHEAP commandline
"""
@with_kw struct _SheapDefaultOptions
    format_in::String = "xyz"
    metadata_file = nothing

    verbose::Bool = false
    quiet::Bool = false
    dim::Int = 0  # Default might change - always specify
    scale::Bool = true
    similarity_threshold::Float64 = 0.0
    cost_threshold::Float64 = 0.0
    use_tsne::Bool = false
    use_umap::Bool = false
    perplexity::Float64 = 15.0
    knn::Int = 15
    uniform_random::Bool = false
    uniform_random_packing::Float64 = 0.5
    compression_factor::Float64 = 10.0
    kl_divergence_loss::Bool = false
    cross_entropy_loss::Bool = true
    hard_sphere::Bool = false
    hard_sphere_core_strength::Float64 = 1.0
    hard_sphere_steps::Int = 500
    hard_sphere_tol::Int = 2.0
    sphere_radius::Float64 = 0.02
end

const SHEAP_OPT_MAP = Dict(
    :format_in => "-read",
    :metadata_file => "-m",
    :verbose => "-v",
    :quiet => "-q",
    :dim => "-dim",
    :scale => "-scale",
    :similarity_threshold => "-st",
    :cost_threshold => "-et",
    :use_tsne => "-tsne",
    :use_umap => "-umap",
    :perplexity => "-p",
    :knn => "-k",
    :uniform_random => "-up",
    :uniform_random_packing => "-f",
    :compression_factor => "-pca",
    :kl_divergence_loss => "-kl",
    :cross_entropy_loss => "-ce",
    :hard_sphere => "-hs",
    :hard_sphere_core_strength => "-cs",
    :hard_sphere_steps => "-gs",
    :hard_sphere_tol => "-grtol",
    :sphere_radius => "-rs",
)

# Check for the consistency of map
let
    optkey = keys(SHEAP_OPT_MAP)
    fields = fieldnames(SheapOptions)
    for name in fields
        @assert name in optkey "SHEAP_OPT_MAP and SheapOptions are not one-to-one mapped"
    end
    for name in optkey
        @assert name in fields "SHEAP_OPT_MAP and SheapOptions are not one-to-one mapped"
    end
end

const SHEAP_OPT_DEFAULT = _SheapDefaultOptions()


function build_cmd(opt::SheapOptions)
    out = String["sheap"]
    for name in fieldnames(SheapOptions)
        val = getproperty(opt, name)
        getproperty(SHEAP_OPT_DEFAULT, name) == val && continue

        if val == false || isnothing(val)
            continue
        end
        if val == true
            push!(out, SHEAP_OPT_MAP[name])
            continue
        end

        push!(out, SHEAP_OPT_MAP[name])
        push!(out, string(val))
    end
    Cmd(out)
end

@with_kw struct SheapMetadata
    label::String = "SHEAP-IN"
    natoms::Int = 1
    form::String = "Al"
    sym::String = "P1"
    volume::Float64 = -0.1
    enthalpy::Float64 = -0.1
    nfound::Int = 1
end

function write_sheap_input!(out, vecs, metadata=repeat([SheapMetadata()], length(vecs)))
    for (m, vec) in zip(metadata, vecs)
        # Write the output string
        write(out, string(length(vec)) * "\n")
        write(out, join(map(string, vec), "\t"))
        write(out, "\n")
        # Write metadata
        metaline = @sprintf "%s\t%d\t%s\t\"%s\"\t%f\t%f\t%d\n" m.label m.natoms m.form m.sym m.volume m.enthalpy m.nfound
        write(out, metaline)
    end
end


"""
    run_sheap(vecs, opt::SheapOptions;metadata=repeat([SheapMetadata()], length(vecs)), show_stderr=true)

Run SHEAP for a iterator of vectors with the given options.
"""
function run_sheap(vecs, opt::SheapOptions; metadata=repeat([SheapMetadata()], length(vecs)), show_stderr=false)

    # Store the output in a temporay buffer
    outbuffer = IOBuffer()
    cmd = build_cmd(opt)
    show_stderr ? stderopt = nothing : stderopt = devnull
    pip = pipeline(cmd, stdout=outbuffer, stderr=stderopt)

    open(pip, "w") do handle
        write_sheap_input!(handle, vecs, metadata)
    end
    seek(outbuffer, 0)
    # Now parse the data
    parse_sheap_output(outbuffer)
end

"""
Prase the output of SHEAP from an IO object
"""
function parse_sheap_output(buffer::IO)
    nitem = 0
    dim = 0
    cost = 0.0
    local coords
    local form
    local labels
    local nforms
    local sym
    local c1
    local c2
    local nfound
    local radius
    offset = 1
    errors = []
    for (i, line) in enumerate(eachline(buffer))
        if offset == i
            try
                nitem = parse(Int, line)
                if !isempty(errors)
                    @warn "Errors in STDOUT: $(join(errors, "\n"))"
                end
            catch
                offset += 1
                push!(errors, line)
            end
            continue
        end

        if i == offset + 1
            tokens = split(line)
            dim = parse(Int, tokens[2])
            cost = parse(Float64, tokens[end])
            coords = Array{Float64,2}(undef, nitem, dim)
            labels = Vector{String}(undef, nitem)
            form = Vector{String}(undef, nitem)
            nforms = Vector{Int}(undef, nitem)
            sym = Vector{String}(undef, nitem)
            c1 = Vector{Float64}(undef, nitem)
            c2 = Vector{Float64}(undef, nitem)
            nfound = Vector{Int}(undef, nitem)
            radius = Vector{Float64}(undef, nitem)
            continue
        end

        # Read the data
        tokens = split(line)
        for d = 1:dim
            coords[i-1-offset, d] = parse(Float64, tokens[d+1])
        end

        labels[i-1-offset] = tokens[2+dim]
        nforms[i-1-offset] = parse(Int, tokens[3+dim])
        form[i-1-offset] = tokens[4+dim]
        sym[i-1-offset] = tokens[5+dim]
        c1[i-1-offset] = parse(Float64, tokens[6+dim])
        c2[i-1-offset] = parse(Float64, tokens[7+dim])
        nfound[i-1-offset] = parse(Int, tokens[8+dim])
        radius[i-1-offset] = parse(Float64, tokens[9+dim])
    end
    (; labels, nforms, form, sym, c1, c2, nfound, radius, coords)
end
end # END module

using .SheapIO