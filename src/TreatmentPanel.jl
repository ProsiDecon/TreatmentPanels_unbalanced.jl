using DataFrames, Dates, Parameters

# Abstract type for dispatch
abstract type TreatmentPanel end

# Types for number of treatment units and periods
# Treatment duration - this is either continuous or discontinuous
abstract type TreatmentDurationType end
struct Continuous <: TreatmentDurationType end
struct Discontinuous <: TreatmentDurationType end

# Treatment timing - relevant only for MultiUnitTreatments, can be simultaneous or staggered
abstract type TreatmentTimingType end
struct Simultaneous{T <: TreatmentDurationType} <: TreatmentTimingType end
struct Staggered{T <: TreatmentDurationType} <: TreatmentTimingType end

# Unit type - either single or multiple treated units
abstract type TreatmentType end
struct SingleUnitTreatment{T <: TreatmentDurationType} <: TreatmentType end
struct MultiUnitTreatment{T <: TreatmentTimingType}  <: TreatmentType end


# BalancedPanel will have an N×T matrix of treatment assigment and outcomes
"""
    BalancedPanel{TreatmentType}

A TreatmentPanel in which all N treatment units are observed for the same T periods.

The object is constructed from a `DataFrame` which contains all outcome and covariate data. The 
constructor requires passing information on which column in the `DataFrame` holds the subject and 
time period identifiers, as well as information on the timing of treatment. Treatments are specified
as a `Pair` of either a `String` or `Symbol` identifying the unit treated, and a `Date` or
`Int` value, identifying the treatment period. Where treatments have an end point, the
treatment period os specified as a `Tuple` of either `Date` or `Int`, indicating the first
and last period of treatment. Where individual units have multiple treatment periods,
these are specified as a `Pair{Union{String, Symbol}, Vector{Union{Int, Date}}}` of treatment
unit identifier and vector of timings. Finally, single continuous, single periodic, and
multiple period treatments can be generalized to multiple treated units.
The following table provides an overview of the types of treatment pattern supported:


|                     |  Only starting point        |   Start and end point                     |   Multiple start & end points                       |
|---------------------|-----------------------------|-------------------------------------------|-----------------------------------------------------|
| **one unit**        |  Pair{String, Date}         |   Pair{String, Tuple{Date, Date}}         |  Pair{String}, Vector{Tuple{Date, Date}}}           |
| **multiple units**  |  Vector{Pair{String, Date}} |   Vector{Pair{String, Tuple{Date, Date}}} |  Vector{Pair{String}, Vector{Tuple{Date, Date}}}}   |
"""
@with_kw struct BalancedPanel{UTType} <: TreatmentPanel where UTType <: TreatmentType
    W::Union{Matrix{Bool}, Matrix{Union{Missing, Bool}}}
    #Y::Union{Matrix{Float64},Array{Float64,3}}
    Y::Array
    df::DataFrame
    id_var::Union{String, Symbol}
    t_var::Union{String, Symbol}
    outcome_var::Union{String, Symbol}
    ts::Vector{T1} where T1 <: Union{Date, Int64}
    is::Vector{T2} where T2 <: Union{Symbol, String, Int64}
end

# Check that ID, time, and outcome variable are provided
function check_id_t_outcome(df, outcome_var, id_var, t_var)
    # Check relevant info has been provided
    !isnothing(outcome_var) || error(ArgumentError(
        "Please specify outcome_var, the name of the column in your dataset holding the "*
        "outcome variable of interest in your dataset."
        ))
    !isnothing(id_var) || error(ArgumentError(
            "Please specify id_var, the name of the column in your data set holding the "*
            "identifier of your units of observation (panel dimension)."
            ))
    !isnothing(t_var) || error(ArgumentError(
            "Please specify t_var, the name of the column in your dataset holding the "*
            "time dimension."
        ))

    # Ensure columns exist in data
    f = in(names(df))
    f(string(id_var)) || throw("Error: ID variable $id_var is not present in the data.")
    f(string(t_var)) || throw("Error: Time variable $t_var is not present in the data.")
    f(string(outcome_var)) || throw("Error: ID variable $outcome_var is not present in the data.")
end

# Functions to get all treatment periods
function treatment_periods(ta::Pair{T1, S1}) where T1 where S1
    [last(ta)]
end

function treatment_periods(ta::Pair{T1, S1}) where T1 where S1 <: Union{Pair{Int, Int}, Pair{Date, Date}}
    collect(last(ta))
end

function treatment_periods(ta::Vector{Pair{T1, S1}}) where T1 where S1
    last.(ta)
end

function treatment_periods(ta::Vector{Pair{T1, S1}}) where T1 where S1 <: Union{Pair{Int, Int}, Pair{Date, Date}}
    unique(reduce(vcat, collect.(last.(ta))))
end

# Functions to construct treatment assignment matrix
function construct_W(ta::Pair{T1, S1}, N, T, is, ts) where T1 where S1
    W = [false for i = 1:N, j = 1:T]
    W[findfirst(==(ta[1]), is), findfirst(==(ta[2]), ts):end] .= true

    return W
end

function construct_W(ta::Pair{T1, S1}, N, T, is, ts) where T1 where S1 <: Union{Pair{Int, Int}, Pair{Date, Date}}
    W = [false for i = 1:N, j = 1:T]
    W[findfirst(==(ta[1]), is), findfirst(==(ta[2][1]), ts):findfirst(==(ta[2][2]), ts)] .= true

    return W
end

function construct_W(tas::Vector{Pair{T1, S1}}, N, T, is, ts) where T1 where S1
    W = [false for i = 1:N, j = 1:T]
    for ta ∈ tas
        W[findfirst(==(ta[1]), is), findfirst(==(ta[2]), ts):end] .= true
    end

    return W
end

function construct_W(tas::Vector{Pair{T1, S1}}, N, T, is, ts) where T1 where S1 <: Union{Pair{Int, Int}, Pair{Date, Date}}
    W = [false for i = 1:N, j = 1:T]
    for ta ∈ tas
        W[findfirst(==(ta[1]), is), findfirst(==(ta[2][1]), ts):findfirst(==(ta[2][2]), ts)] .= true
    end

    return W
end


## DP: Version with covariates 
# Constructor for single continuous treatment - returns BalancedPanel{SingleUnitTreatment{Continuous}}
function BalancedPanel_maker(df::DataFrame, 
                        treatment_assignment::Pair{T1, T2};
                        id_var = nothing, 
                        t_var = nothing, 
                        outcome_var = nothing, 
                        covariates::Union{Nothing, Vector{Symbol}, Symbol} = nothing,
                        sort_inplace = false) where T1 where T2 <: Union{Date, Int}

    if typeof(covariates) == Symbol
        covariates = [covariates]
    end
                        
    # Get all units and time periods
    is = sort(unique(df[!, id_var])); i_set = Set(is)
    ts = sort(unique(df[!, t_var])); t_set = Set(ts)

    # Dimensions
    N = length(is)
    T = length(ts)

    # Get all treatment units and treatment periods
    treated_i = first(treatment_assignment)
    treated_t = last(treatment_assignment)

    # Sanity checks
    check_id_t_outcome(df, outcome_var, id_var, t_var)
    in(treated_i, i_set) || throw("Error: Treatment unit $treated_i is not in the list of unit identifiers $id_var")
    in(treated_t, t_set) || throw("Error: Treatment period $treated_t is not in the list of time identifiers $t_var")
    
    # Sort data if necessary, in place if required
    df = ifelse(issorted(df, [id_var, t_var]), df, 
                                               ifelse(sort_inplace, sort!(df, [id_var, t_var]), 
                                                                    sort(df, [id_var, t_var])))

    # Treatment matrix
    W = construct_W(treatment_assignment, N, T, is, ts)
    #W = TreatmentPanels.construct_W(treatment_assignment, N, T, is, ts)

    # Outcome matrix
    if !isnothing(covariates)
        Y = zeros(eltype(df[!, outcome_var]), (size(W)..., length(covariates)+1))
        for (row, i) ∈ enumerate(is), (col, t) ∈ enumerate(ts)
            Y[row, col,1] = only(df[(df[!, id_var] .== i) .& (df[!, t_var] .== t), outcome_var])
        end
        # add covariates in 3rd dimension
        for covariate in eachindex(covariates)
            for (row, i) ∈ enumerate(is), (col, t) ∈ enumerate(ts)
                Y[row, col,1+covariate] = only(df[(df[!, id_var] .== i) .& (df[!, t_var] .== t), covariates[covariate]])
            end
        end
    else
        Y = zeros(eltype(df[!, outcome_var]), size(W))
        for (row, i) ∈ enumerate(is), (col, t) ∈ enumerate(ts)
            Y[row, col] = only(df[(df[!, id_var] .== i) .& (df[!, t_var] .== t), outcome_var])
        end
    end

    return BalancedPanel{SingleUnitTreatment{Continuous}}(W, Y, df, id_var, t_var, outcome_var, ts, is)  
end


# Getter functions
"""
    treated_ids(x <: BalancedPanel)

    Returns the indices of treated units in the panel, so that Y[treated_ids(x), :] returns a
    (Nₜᵣ×T) matrix of outcomes for treated units in all periods.
"""
function treated_ids(x::BalancedPanel{SingleUnitTreatment{T}}) where T
    for i ∈ 1:size(x.Y, 1)
        for t ∈ 1:size(x.Y, 2)
            if x.W[i, t]
                return i
            end
        end
    end
end

function treated_ids(x::BalancedPanel{MultiUnitTreatment{T}}) where T
    any.(eachrow(x.W))
end


function treated_ids(x::BalancedPanelQ{SingleUnitTreatment{T}}) where T
    for i ∈ 1:size(x.Y, 1)
        for t ∈ 1:size(x.Y, 2)
            if x.W[i, t]
                return i
            end
        end
    end
end

function treated_ids(x::BalancedPanelQ{MultiUnitTreatment{T}}) where T
    any.(eachrow(x.W))
end


"""
    treated_labels(x <: BalancedPanel)

    Returns the labels of treated units as given by the `id_var` column in the underlying data set.
"""
function treated_labels(x::BalancedPanel{SingleUnitTreatment{T}}) where T
    x.is[treated_ids(x)]
end

function treated_labels(x::BalancedPanel{MultiUnitTreatment{Simultaneous{Continuous}}})
    x.is[treated_ids(x)]
end

function treated_labels(x::BalancedPanelQ{SingleUnitTreatment{T}}) where T
    x.is[treated_ids(x)]
end

function treated_labels(x::BalancedPanelQ{MultiUnitTreatment{Simultaneous{Continuous}}})
    x.is[treated_ids(x)]
end

"""
    first_treated_period_ids(x <: BalancedPanel)

    Returns the indices of the first treated period for each treated units, that is, a Vector{Int}
    of length Nₜᵣ, where each element is the index of the first 1 in the row of treatment matrix W
    corresonding to the treatment unit. 

    For a single treated unit, returns only the index of the first treated period. For multiple 
"""
function first_treated_period_ids(x::BalancedPanel{SingleUnitTreatment{T}}) where T
    findfirst(vec(x.W[treated_ids(x), :]))
end

function first_treated_period_ids(x::BalancedPanel{MultiUnitTreatment{T}}) where T
    [x for x ∈ findfirst.(eachrow(x.W)) if !isnothing(x)]
end

"""
    first_treated_period_labels(x <: BalancedPanel)

    Returns the labels of the first treated period for each treated units, that is, a Vector{T}
    of length Nₜᵣ, where T is the eltype of the `t_var` column in the underlying data.
"""
function first_treated_period_labels(x::BalancedPanel{SingleUnitTreatment{T}}) where T
    x.ts[first_treated_period_ids(x)]
end

"""
    length_T₀(x <: BalancedPanel)

    Returns the number of pre-treatment periods for each treated unit.
"""
function length_T₀(x::BalancedPanel{SingleUnitTreatment{Continuous}})
    first_treated_period_ids(x) - 1
end

function length_T₀(x::BalancedPanel{MultiUnitTreatment{Simultaneous{Continuous}}})
    first_treated_period_ids(x) .- 1
end

function length_T₀(x::BalancedPanelQ{SingleUnitTreatment{Continuous}})
    first_treated_period_ids(x) - 1
end

function length_T₀(x::BalancedPanelQ{MultiUnitTreatment{Simultaneous{Continuous}}})
    first_treated_period_ids(x) .- 1
end


"""
    length_T₁(x <: BalancedPanel)

    Returns the number of treatment periods for each treated unit. 
"""
function length_T₁(x::BalancedPanel{SingleUnitTreatment{Continuous}})
    size(x.Y, 2) .- first_treated_period_ids(x) + 1
end

function length_T₁(x::BalancedPanel{MultiUnitTreatment{Simultaneous{Continuous}}})
    size(x.Y[treated_ids(x), :], 2) .- first_treated_period_ids(x) .+ 1
end

""" 
    get_y₁₀(x <: BalancedPanel)

    Returns the pre-treatment outcomes for the treated unit(s). For SingleUnitTreatment designs,
    this is a vector of length T₀, while for MultiUnitTreatment designs, it is a (Nₜᵣ×T₀) matrix
"""
function get_y₁₀(x::BalancedPanel{SingleUnitTreatment{Continuous}})
    x.Y[treated_ids(x), 1:first_treated_period_ids(x)-1]
end


""" 
    get_y₁₁(x <: BalancedPanel)

    Returns the post-treatment outcomes for the treated unit(s). For SingleUnitTreatment designs,
    this is a vector of length T₁, while for MultiUnitTreatment designs, it is a (Nₜᵣ×T₁) matrix
"""
function get_y₁₁(x::BalancedPanel{SingleUnitTreatment{Continuous}})
    x.Y[x.W]
end

""" sc
    get_y₀₀(x <: BalancedPanel)

    Returns the pre-treatment outcomes for the untreated units, an (Nₖₒ×T₀) matrix
"""
function get_y₀₀(x::BalancedPanel{SingleUnitTreatment{Continuous}})
    x.Y[Not(treated_ids(x)), 1:first_treated_period_ids(x)-1]
end

""" 
    get_y₀₀(x <: BalancedPanel)

    Returns the post-treatment outcomes for the untreated units, an (Nₖₒ×T₁) matrix
"""
function get_y₀₁(x::BalancedPanel{SingleUnitTreatment{Continuous}})
    x.Y[Not(treated_ids(x)), first_treated_period_ids(x):end]
end

""" 
    decompose_y(x <: BalancedPanel)

    Decomposes the outcome matrix Y into four elements:
    
    * Pre-treatment outcomes for treated units (y₁₀)
    * Post-treatment outcomes for treated units (y₁₁)
    * Pre-treatment outcomes for control units (y₀₀)
    * Post-treatment outcomes for treated units (y₀₁)

    and returns a tuple (y₁₀, y₁₁, y₀₀, y₀₁)
"""
function decompose_y(x)
    get_y₁₀(x), get_y₁₁(x), get_y₀₀(x), get_y₀₁(x)
end

####################################################################################################

# Constructor for single start/end treatment - returns BalancedPanel{SingleUnitTreatment{Discontinuous}}
function BalancedPanel_maker(df::DataFrame, treatment_assignment::Pair{T1, T2};
    id_var = nothing, t_var = nothing, outcome_var = nothing, 
    sort_inplace = false) where T1 where T2 <: Union{Pair{Date, Date}, Pair{Int, Int}}

    # Get all units and time periods
    is = sort(unique(df[!, id_var])); i_set = Set(is)
    ts = sort(unique(df[!, t_var])); t_set = Set(ts)

    # Dimensions
    N = length(is)
    T = length(ts)

    # Get all treatment units and treatment periods
    treated_i = first(treatment_assignment)
    treatment_start_end = last(treatment_assignment)

    # Sanity checks
    check_id_t_outcome(df, outcome_var, id_var, t_var)
    in(treated_i, i_set) || throw("Error: Treatment unit $treated_i is not in the list of unit identifiers $id_var")
    in(treatment_start_end[1], t_set) || throw("Error: Treatment start point $(treatment_start_end[1]) is not in the list of time identifiers $t_var")
    in(treatment_start_end[2], t_set) || throw("Error: Treatment end point $(treatment_start_end[2]) is not in the list of time identifiers $t_var")
    
    # Sort data if necessary, in place if required
    df = ifelse(issorted(df, [id_var, t_var]), df, 
                                               ifelse(sort_inplace, sort!(df, [id_var, t_var]), 
                                                                    sort(df, [id_var, t_var])))

    # Treatment matrix
    W = construct_W(treatment_assignment, N, T, is, ts)

    # Outcome matrix
    Y = zeros(size(W))
    for (row, i) ∈ enumerate(is), (col, t) ∈ enumerate(ts)
        try
            Y[row, col] = only(df[(df[!, id_var] .== i) .& (df[!, t_var] .== t), outcome_var])
        catch ArgumentError
            throw("$(nrow(df[(df[!, id_var] .== i) .& (df[!, t_var] .== t), :])) outcomes present in the data for unit $i in period $t")
        end
    end

    BalancedPanel{SingleUnitTreatment{Discontinuous}}(W, Y, df, id_var, t_var, outcome_var, ts, is)  
end


# DP Fallback method with covariates - if the length of treatment assignment is one use single treatment method above
function BalancedPanel_maker(df::DataFrame, treatment_assignment::Union{Pair{T1, Pair{T3, T4}}, Vector{Pair{Int64, Int64}}};
                                id_var = nothing, 
                                t_var = nothing, 
                                outcome_var = nothing, 
                                covariates::Union{Nothing, Vector{Symbol}, Symbol} = nothing,
                                sort_inplace = false)  where {T1, T3, T4} 

    if length(treatment_assignment) == 1
        return BalancedPanel_maker(df, only(treatment_assignment); id_var = id_var, t_var = t_var,
                                outcome_var = outcome_var, covariates = covariates, sort_inplace = sort_inplace)
    end

    if typeof(covariates) == Symbol
        covariates = [covariates]
    end

    # Get all units and time periods
    is = sort(unique(df[!, id_var])); i_set = Set(is)
    ts = sort(unique(df[!, t_var])); t_set = Set(ts)

    # Get all treatment units and treatment periods
    treated_is = first.(treatment_assignment)
    treated_is = typeof(treated_is) <: AbstractArray ? treated_is : [treated_is]
    treated_ts = treatment_periods(treatment_assignment)
    #treated_ts = TreatmentPanels.treatment_periods(treatment_assignment)

    # Dimensions
    N = length(is)
    T = length(ts)

    ### SANITY CHECKS ###
    #TreatmentPanels.check_id_t_outcome(df, outcome_var, id_var, t_var)
    check_id_t_outcome(df, outcome_var, id_var, t_var)
    for ti ∈ treated_is
        in(ti, i_set) || throw("Error: Treatment unit $ti is not in the list of unit identifiers $id_var")
    end

    for tt ∈ treated_ts
        in(tt, t_set) || throw("Error: Treatment period $tt is not in the list of time identifiers $t_var")
    end
    
    # Sort data if necessary, in place if required
    df = ifelse(issorted(df, [id_var, t_var]), df, 
                                               ifelse(sort_inplace, sort!(df, [id_var, t_var]), 
                                                                    sort(df, [id_var, t_var])))

    # Treatment matrix
    W = construct_W(treatment_assignment, N, T, is, ts)
    #W = TreatmentPanels.construct_W(treatment_assignment, N, T, is, ts)
    
     # Outcome matrix
     if !isnothing(covariates)
        Y = zeros(eltype(df[!, outcome_var]), (size(W)..., length(covariates)+1))
        for (row, i) ∈ enumerate(is), (col, t) ∈ enumerate(ts)
            Y[row, col,1] = only(df[(df[!, id_var] .== i) .& (df[!, t_var] .== t), outcome_var])
        end
        # add covariates in 3rd dimension
        for covariate in eachindex(covariates)
            for (row, i) ∈ enumerate(is), (col, t) ∈ enumerate(ts)
                Y[row, col,1+covariate] = only(df[(df[!, id_var] .== i) .& (df[!, t_var] .== t), covariates[covariate]])
            end
        end
    else
        Y = zeros(eltype(df[!, outcome_var]), size(W))
        for (row, i) ∈ enumerate(is), (col, t) ∈ enumerate(ts)
            Y[row, col] = only(df[(df[!, id_var] .== i) .& (df[!, t_var] .== t), outcome_var])
        end
    end

    # Determine TreatmentType and TreatmentDurationType
    uttype = if all(==(treatment_assignment[1][2]), last.(treatment_assignment))
        Simultaneous
    else
        Staggered
    end

    tdtype = if typeof(treatment_assignment) <: Pair
        if typeof(treatment_assignment[2]) <: Pair
            Discontinuous
        else
            Continuous
        end
    else
        if typeof(treatment_assignment[1][2]) <: Pair
            Discontinuous
        else
            Continuous
        end
    end

    BalancedPanel{MultiUnitTreatment{uttype{tdtype}}}(W, Y, df, id_var, t_var, outcome_var, ts, is)  
end



###################################################################################################
### add Balanced Panel structs with baseline weights Q (either a dyadic tensor object or an N × T matrix)
@with_kw struct BalancedPanelQ{UTType} <: TreatmentPanel where UTType <: TreatmentType
    W::Union{Matrix{Bool}, Matrix{Union{Missing, Bool}}}
    #Y::Union{Matrix{Float64},Array{Float64,3}}
    Y::Array
    Q::Array                # the baseline weights 
    df::DataFrame
    id_var::Union{String, Symbol}
    t_var::Union{String, Symbol}
    outcome_var::Union{String, Symbol}
    ts::Vector{T1} where T1 <: Union{Date, Int64}
    is::Vector{T2} where T2 <: Union{Symbol, String, Int64}
end

# Constructor for single continuous treatment - returns BalancedPanel{SingleUnitTreatment{Continuous}}
function BalancedPanelQ_maker(df::DataFrame, 
                        treatment_assignment::Pair{<:Any, <:Union{Date, Int}},
                        df_q::Union{Nothing, DataFrame} = nothing;
                        id_var = nothing, 
                        t_var = nothing, 
                        outcome_var = nothing, 
                        q_var = nothing,                # the variable name if the weights in either df or df_q
                        q_id_var = nothing,             # if q_dyadic, the variable name of the (treated) partner of the current weight
                        q_dyadic::Bool = true,
                        covariates::Union{Nothing, Vector{Symbol}, Symbol} = nothing,
                        sort_inplace = false) #where T1 where T2 <: Union{Date, Int}

    if typeof(covariates) == Symbol
        covariates = [covariates]
    end
                        
    # Get all units and time periods
    is = sort(unique(df[!, id_var])); i_set = Set(is)
    ts = sort(unique(df[!, t_var])); t_set = Set(ts)

    # Dimensions
    N = length(is)
    T = length(ts)

    # Get all treatment units and treatment periods
    treated_i = first(treatment_assignment)
    treated_t = last(treatment_assignment)

    # Sanity checks
    check_id_t_outcome(df, outcome_var, id_var, t_var)
    in(treated_i, i_set) || throw("Error: Treatment unit $treated_i is not in the list of unit identifiers $id_var")
    in(treated_t, t_set) || throw("Error: Treatment period $treated_t is not in the list of time identifiers $t_var")
    
    # Sort data if necessary, in place if required
    df = ifelse(issorted(df, [id_var, t_var]), df, 
                                               ifelse(sort_inplace, sort!(df, [id_var, t_var]), 
                                                                    sort(df, [id_var, t_var])))

    # Treatment matrix
    W = construct_W(treatment_assignment, N, T, is, ts)
    #W = TreatmentPanels.construct_W(treatment_assignment, N, T, is, ts)

    # Outcome matrix
    if !isnothing(covariates)
        Y = zeros(eltype(df[!, outcome_var]), (size(W)..., length(covariates)+1))
        for (row, i) ∈ enumerate(is), (col, t) ∈ enumerate(ts)
            Y[row, col,1] = only(df[(df[!, id_var] .== i) .& (df[!, t_var] .== t), outcome_var])
        end
        # add covariates in 3rd dimension
        for covariate in eachindex(covariates)
            for (row, i) ∈ enumerate(is), (col, t) ∈ enumerate(ts)
                Y[row, col,1+covariate] = only(df[(df[!, id_var] .== i) .& (df[!, t_var] .== t), covariates[covariate]])
            end
        end
    else
        Y = zeros(eltype(df[!, outcome_var]), size(W))
        for (row, i) ∈ enumerate(is), (col, t) ∈ enumerate(ts)
            Y[row, col] = only(df[(df[!, id_var] .== i) .& (df[!, t_var] .== t), outcome_var])
        end
    end

    # make the object for baseline weights
    if q_dyadic 
        !isnothing(df_q) || error(ArgumentError("Dyadic weights chosen but no dataframe df_q for weights provided."))
        !isnothing(q_var) || error(ArgumentError(
            "Please specify q_var, the name of the column in your dataset holding the weight variable in your weight dataset df_q."))
        !isnothing(q_id_var) || error(ArgumentError(
            "Please specify q_id_var, the name of the column in your dataset holding the dyadic partner for the weight in your weight dataset df_q. Otherwise specify q_dyadic = false"))

        Q = zeros(eltype(df_q[!, q_var]), (size(W)..., N))

        q_lookup = Dict{Tuple{eltype(is), eltype(is), eltype(ts)}, eltype(df_q[!, q_var])}()

        for row in eachrow(df_q)
            key = (row[id_var], row[q_id_var], row[t_var])
            q_lookup[key] = row[q_var]  # Overwrites duplicates, change if needed
        end

        # Now loop and fill Q efficiently
        for (row, i) ∈ enumerate(is), (col, t) ∈ enumerate(ts), (slice, j) ∈ enumerate(is)
            Q[row, col, slice] = get(q_lookup, (i, j, t), missing)
        end

        #for (row, i) ∈ enumerate(is), (col, t) ∈ enumerate(ts), (slice, j) ∈ enumerate(is)
        #    Q[row, col, slice] = only(df_q[(df_q[!, id_var] .== i) .& (df_q[!, q_id_var] .== j) .& (df_q[!, t_var] .== t), q_var])
        #end
    else
        if isnothing(df_q) && !isnothing(q_var)
            df_q = select(df, [id_var, t_var, q_var])
        elseif isnothing(df_q) && isnothing(q_var)
            error(ArgumentError("Neither dataframe for weights df_q provided nor variable q_var for weights specified."))
        end

        Q = zeros(eltype(df[!, q_var]), (size(W)))
        for (row, i) ∈ enumerate(is), (col, t) ∈ enumerate(ts)
            Q[row, col] = only(df_q[(df_q[!, id_var] .== i) .& (df_q[!, t_var] .== t), q_var])
        end
    end

    return BalancedPanelQ{SingleUnitTreatment{Continuous}}(W, Y, Q, df, id_var, t_var, outcome_var, ts, is)  
end




# Fallback method - if the length of treatment assignment is one use single treatment method above
function BalancedPanelQ_maker(df::DataFrame, treatment_assignment,
                                df_q::Union{Nothing, DataFrame} = nothing;
                                id_var = nothing, 
                                t_var = nothing, 
                                outcome_var = nothing, 
                                q_var = nothing,
                                q_id_var = nothing,
                                q_dyadic::Bool = true,
                                covariates::Union{Nothing, Vector{Symbol}, Symbol} = nothing,
                                sort_inplace = false)

    if length(treatment_assignment) == 1
        return BalancedPanel_maker(df, only(treatment_assignment); id_var = id_var, t_var = t_var,
                                outcome_var = outcome_var, covariates = covariates, sort_inplace = sort_inplace)
    end

    if typeof(covariates) == Symbol
        covariates = [covariates]
    end

    # Get all units and time periods
    is = sort(unique(df[!, id_var])); i_set = Set(is)
    ts = sort(unique(df[!, t_var])); t_set = Set(ts)

    # Get all treatment units and treatment periods
    treated_is = first.(treatment_assignment)
    treated_is = typeof(treated_is) <: AbstractArray ? treated_is : [treated_is]
    treated_ts = treatment_periods(treatment_assignment)
    #treated_ts = TreatmentPanels.treatment_periods(treatment_assignment)

    # Dimensions
    N = length(is)
    T = length(ts)

    ### SANITY CHECKS ###
    #TreatmentPanels.check_id_t_outcome(df, outcome_var, id_var, t_var)
    check_id_t_outcome(df, outcome_var, id_var, t_var)
    for ti ∈ treated_is
        in(ti, i_set) || throw("Error: Treatment unit $ti is not in the list of unit identifiers $id_var")
    end

    for tt ∈ treated_ts
        in(tt, t_set) || throw("Error: Treatment period $tt is not in the list of time identifiers $t_var")
    end
    
    # Sort data if necessary, in place if required
    df = ifelse(issorted(df, [id_var, t_var]), df, 
                                               ifelse(sort_inplace, sort!(df, [id_var, t_var]), 
                                                                    sort(df, [id_var, t_var])))

    # Treatment matrix
    W = construct_W(treatment_assignment, N, T, is, ts)
    #W = TreatmentPanels.construct_W(treatment_assignment, N, T, is, ts)
    
    # Outcome matrix
    if !isnothing(covariates)
        Y = zeros(eltype(df[!, outcome_var]), (size(W)..., length(covariates)+1))
        for (row, i) ∈ enumerate(is), (col, t) ∈ enumerate(ts)
            Y[row, col,1] = only(df[(df[!, id_var] .== i) .& (df[!, t_var] .== t), outcome_var])
        end
        # add covariates in 3rd dimension
        for covariate in eachindex(covariates)
            for (row, i) ∈ enumerate(is), (col, t) ∈ enumerate(ts)
                Y[row, col,1+covariate] = only(df[(df[!, id_var] .== i) .& (df[!, t_var] .== t), covariates[covariate]])
            end
        end
    else
        Y = zeros(eltype(df[!, outcome_var]), size(W))
        for (row, i) ∈ enumerate(is), (col, t) ∈ enumerate(ts)
            Y[row, col] = only(df[(df[!, id_var] .== i) .& (df[!, t_var] .== t), outcome_var])
        end
    end

    # make the object for baseline weights
    if q_dyadic 
        !isnothing(df_q) || error(ArgumentError("Dyadic weights chosen but no dataframe df_q for weights provided."))
        !isnothing(q_var) || error(ArgumentError(
            "Please specify q_var, the name of the column in your dataset holding the weight variable in your weight dataset df_q."))
        !isnothing(q_id_var) || error(ArgumentError(
            "Please specify q_id_var, the name of the column in your dataset holding the dyadic partner for the weight in your weight dataset df_q. Otherwise specify q_dyadic = false"))

        Q = zeros(eltype(df_q[!, q_var]), (size(W)..., N))

        q_lookup = Dict{Tuple{eltype(is), eltype(is), eltype(ts)}, eltype(df_q[!, q_var])}()

        for row in eachrow(df_q)
            key = (row[id_var], row[q_id_var], row[t_var])
            q_lookup[key] = row[q_var]  # Overwrites duplicates, change if needed
        end

        # Now loop and fill Q efficiently
        for (row, i) ∈ enumerate(is), (col, t) ∈ enumerate(ts), (slice, j) ∈ enumerate(is)
            Q[row, col, slice] = get(q_lookup, (i, j, t), missing)
        end

        #for (row, i) ∈ enumerate(is), (col, t) ∈ enumerate(ts), (slice, j) ∈ enumerate(is)
        #    Q[row, col, slice] = only(df_q[(df_q[!, id_var] .== i) .& (df_q[!, q_id_var] .== j) .& (df_q[!, t_var] .== t), q_var])
        #end
    else
        if isnothing(df_q) && !isnothing(q_var)
            df_q = select(df, [id_var, t_var, q_var])
        elseif isnothing(df_q) && isnothing(q_var)
            error(ArgumentError("Neither dataframe for weights df_q provided nor variable q_var for weights specified."))
        end

        Q = zeros(eltype(df[!, q_var]), (size(W)))
        for (row, i) ∈ enumerate(is), (col, t) ∈ enumerate(ts)
            Q[row, col] = only(df_q[(df_q[!, id_var] .== i) .& (df_q[!, t_var] .== t), q_var])
        end
    end

    # Determine TreatmentType and TreatmentDurationType
    uttype = if all(==(treatment_assignment[1][2]), last.(treatment_assignment))
        Simultaneous
    else
        Staggered
    end

    tdtype = if typeof(treatment_assignment) <: Pair
        if typeof(treatment_assignment[2]) <: Pair
            Discontinuous
        else
            Continuous
        end
    else
        if typeof(treatment_assignment[1][2]) <: Pair
            Discontinuous
        else
            Continuous
        end
    end

    BalancedPanelQ{MultiUnitTreatment{uttype{tdtype}}}(W, Y, Q, df, id_var, t_var, outcome_var, ts, is)  
end


## UnblancedPanel - N observations but not all of them for T periods

#!# Not yet implemented