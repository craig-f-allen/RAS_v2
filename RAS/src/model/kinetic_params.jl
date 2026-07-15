using XLSX
using Base.Iterators: drop

# Warning to yee who enter here: some overloading trickery afoot in order to solve the problem of making one clean simulation using two different ways of dynamically calculating parameters. There may be simpler ways than this but this is type stable and clean.

abstract type AbstractKineticParams end

# ----------------------------------------------------------
# Struct 1: Wild Type Layout (kT is virtual)
# ----------------------------------------------------------
mutable struct WTKineticParams <: AbstractKineticParams
    kint::Float64
    kdGDP::Float64
    kdGTP::Float64
    kaGDP::Float64
    kaGTP::Float64    # Stored parameter
    kcat::Float64
    KM::Float64
    kGDP::Float64
    KMGDP::Float64
    KMGTP::Float64
    # kGTP is virtual
    KD::Float64
    kaEff::Float64
    # kdEff is virtual
end

function Base.getproperty(obj::WTKineticParams, sym::Symbol)
    val::Float64 = if sym === :kGTP
        # kGTP = kGDP * KMGTP * ((kaGDP * kdGTP) / (kdGDP * kaGTP)) / KMGDP
        return getfield(obj, :kGDP) * getfield(obj, :KMGTP) * 
               ((getfield(obj, :kaGDP) * getfield(obj, :kdGTP)) / 
                (getfield(obj, :kdGDP) * getfield(obj, :kaGTP))) / getfield(obj, :KMGDP)
                
    elseif sym === :kdEff
        # kdEff = kaEff * kGDP
        getfield(obj, :kaEff) * getfield(obj, :kGDP)
        
    else
        getfield(obj, sym)
    end
    return val
end

# Allows autocomplete and introspection tools to see the virtual properties
Base.propertynames(::WTKineticParams) = (fieldnames(WTKineticParams)..., :kGTP, :kdEff)

# ------------------------------------------------------------------------------
# Struct 2: Mutant Layout (kaGTP is virtual/calculated)
# ------------------------------------------------------------------------------
mutable struct MutantKineticParams <: AbstractKineticParams
    kint::Float64
    kdGDP::Float64
    kdGTP::Float64
    kaGDP::Float64
    # kaGTP is virtual
    kcat::Float64
    KM::Float64
    kGDP::Float64
    KMGDP::Float64
    KMGTP::Float64
    kGTP::Float64     # Stored parameter
    KD::Float64
    kaEff::Float64
    # kdEff is virtual
end

function Base.getproperty(obj::MutantKineticParams, sym::Symbol)
    val::Float64 = if sym === :kaGTP
        # kaGTP = kGDP * KMGTP * ((kaGDP * kdGTP) / (kdGDP * kGTP)) / KMGDP
        return getfield(obj, :kGDP) * getfield(obj, :KMGTP) * 
               ((getfield(obj, :kaGDP) * getfield(obj, :kdGTP)) / 
                (getfield(obj, :kdGDP) * getfield(obj, :kGTP))) / getfield(obj, :KMGDP)
                
    elseif sym === :kdEff
        # kdEff = kaEff * kGDP
        getfield(obj, :kaEff) * getfield(obj, :kGDP)
        
    else
        getfield(obj, sym)
    end
    return val
end

# ------------------------------------------------------------------------------
# Initialization & Data Loading
# ------------------------------------------------------------------------------
const vol_scale = 250

# Initialized exactly in your specified order (excluding the virtual property :kGTP)
WT = WTKineticParams(
    3.5e-4,             # 1.  kint
    1.1e-4,             # 2.  kdGDP
    2.5e-4,             # 3.  kdGTP
    2.3e6,              # 4.  kaGDP
    2.2e6,              # 5.  kaGTP
    5.4,                # 6.  kcat
    0.23e-6/vol_scale,  # 7.  KM
    3.9,                # 8.  kGDP
    3.86e-4/vol_scale,  # 9.  KMGDP
    3e-4/vol_scale,     # 10. KMGTP
    80e-9,              # 12. KD
    4.5e7               # 13. kaEff
)

function get_mutant_params(file_path::String)
    mutants = Dict{Symbol, MutantKineticParams}()
    XLSX.openxlsx(file_path) do f
        sheet = f["Sheet1"]
        for r in drop(XLSX.eachrow(sheet), 1)
            mutants[Symbol(r[1])] = MutantKineticParams(
                convert(Float64, r[2]) * WT.kint,     # 1.  kint
                convert(Float64, r[3]) * WT.kdGDP,    # 2.  kdGDP
                convert(Float64, r[4]) * WT.kdGTP,    # 3.  kdGTP
                convert(Float64, r[5]) * WT.kaGDP,    # 4.  kaGDP
                convert(Float64, r[7]) * WT.kcat,     # 6.  kcat
                convert(Float64, r[8]) * WT.KM,       # 7.  KM
                convert(Float64, r[9]) * WT.kGDP,     # 8.  kGDP
                convert(Float64, r[10]) * WT.KMGDP,   # 9.  KMGDP
                convert(Float64, r[11]) * WT.KMGTP,   # 10. KMGTP
                convert(Float64, r[12]) * WT.kGTP,    # 11. kGTP
                convert(Float64, r[13]) * WT.KD,      # 12. KD
                convert(Float64, r[14]) * WT.kaEff    # 13. kaEff
            )
        end
        return mutants
    end
end
