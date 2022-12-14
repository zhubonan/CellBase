#=
Elemental data
=#
import Base

struct SMACTDataEntry
    symbol::Symbol
    name::String
    Z::Int
    mass::Float64
    r_cov::Float64
    e_affinity::Float64
    p_eig::Float64
    s_eig::Float64
    abundance::Float64
    el_neg::Float64
    ion_pot::Float64
    dipol::Float64
end

function _parse_non_as_nan(x) 
    if x == "None"
        return NaN
    else
        return parse(Float64, x)
    end
end
struct SMACTData
    entries::Vector{SMACTDataEntry}
    by_symbol::Dict{Symbol, SMACTDataEntry}
    by_Z::Dict{Int, SMACTDataEntry}
    by_name::Dict{String, SMACTDataEntry}
end


function SMACTData(entries::Vector)
    by_symbol = Dict{Symbol, SMACTDataEntry}(
            (entry.symbol, entry) for entry in entries
    )
    by_Z = Dict{Int, SMACTDataEntry}(
            (entry.Z, entry) for entry in entries
    )
    by_name = Dict{String, SMACTDataEntry}(
            (entry.name, entry) for entry in entries
    )
    SMACTData(
        entries,
        by_symbol,
        by_Z,
        by_name
    )
end

function _load_smact_data()
    entries = open(joinpath(splitdir(@__FILE__)[1], "element_data.txt")) do fh
        output = SMACTDataEntry[]
        for line in eachline(fh)
            startswith(line, "#") && continue
            tokens = split(line)
            entry = SMACTDataEntry(
                Symbol(tokens[1]),
                tokens[2],
                parse(Int, tokens[3]),
                _parse_non_as_nan(tokens[4]),
                _parse_non_as_nan(tokens[5]),
                _parse_non_as_nan(tokens[6]),
                _parse_non_as_nan(tokens[7]),
                _parse_non_as_nan(tokens[8]),
                _parse_non_as_nan(tokens[9]),
                _parse_non_as_nan(tokens[10]),
                _parse_non_as_nan(tokens[11]),
                _parse_non_as_nan(tokens[12]),
            )
            push!(output, entry)
        end
        output
    end
    SMACTData(entries)
end

const smact_data = _load_smact_data()

Base.getindex(x::SMACTData, key::Int) = x.by_Z[key]
Base.getindex(x::SMACTData, key::Symbol) = x.by_symbol[key]
Base.getindex(x::SMACTData, key::String) = x.by_name[uppercasefirst(key)]