using TreatmentPanels, DataFrames, Parameters, Dates

import TreatmentPanels: treated_ids, treated_labels, first_treated_period_ids, first_treated_period_labels
import TreatmentPanels: length_T₀, length_T₁, get_y₁₀, construct_W, treatment_periods, check_id_t_outcome


###################################################################################################
### add Balanced Panel structs with covariates

@with_kw struct BalancedPanelCov{UTType} <: TreatmentPanel where UTType <: TreatmentType
    W::Union{Matrix{Bool}, Matrix{Union{Missing, Bool}}}
    Y::Array
    df::DataFrame
    id_var::Union{String, Symbol}
    t_var::Union{String, Symbol}
    outcome_var::Union{String, Symbol}
    ts::Vector{T1} where T1 <: Union{Date, Int64}
    is::Vector{T2} where T2 <: Union{Symbol, String, Int64}
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

# function to elicit and promote type of outcome and covariate variables for construction of Y array in BalancedPanelCov
function panel_y_eltype(df::DataFrame, outcome_var, covariates)
    y_vars = [outcome_var; covariates]
    y_eltypes = eltype.(eachcol(df[!, y_vars]))
    all(T -> T <: Real, y_eltypes) ||
        throw(ArgumentError("Outcome and covariates must be numeric."))

    promote_type(y_eltypes...)
end

### multiple dispatch for treated_ids 
function treated_ids(x::BalancedPanelCov{SingleUnitTreatment{T}}) where T
    for i ∈ 1:size(x.Y, 1)
        for t ∈ 1:size(x.Y, 2)
            if x.W[i, t]
                return i
            end
        end
    end
end

function treated_ids(x::BalancedPanelCov{MultiUnitTreatment{T}}) where T
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

# multiple dispatch for treated_labels 

function treated_labels(x::BalancedPanelCov{SingleUnitTreatment{T}}) where T
    x.is[treated_ids(x)]
end

function treated_labels(x::BalancedPanelCov{MultiUnitTreatment{Simultaneous{Continuous}}})
    x.is[treated_ids(x)]
end

function treated_labels(x::BalancedPanelQ{SingleUnitTreatment{T}}) where T
    x.is[treated_ids(x)]
end

function treated_labels(x::BalancedPanelQ{MultiUnitTreatment{Simultaneous{Continuous}}})
    x.is[treated_ids(x)]
end

### multiple dispatch for first treated period ids

function first_treated_period_ids(x::BalancedPanelCov{SingleUnitTreatment{T}}) where T
    findfirst(vec(x.W[treated_ids(x), :]))
end

function first_treated_period_ids(x::BalancedPanelCov{MultiUnitTreatment{T}}) where T
    [x for x ∈ findfirst.(eachrow(x.W)) if !isnothing(x)]
end

function first_treated_period_ids(x::BalancedPanelQ{SingleUnitTreatment{T}}) where T
    findfirst(vec(x.W[treated_ids(x), :]))
end

function first_treated_period_ids(x::BalancedPanelQ{MultiUnitTreatment{T}}) where T
    [x for x ∈ findfirst.(eachrow(x.W)) if !isnothing(x)]
end

# multiple dispatch for first treated period labels
function first_treated_period_labels(x::BalancedPanelCov{SingleUnitTreatment{T}}) where T
    x.ts[first_treated_period_ids(x)]
end

function first_treated_period_labels(x::BalancedPanelQ{SingleUnitTreatment{T}}) where T
    x.ts[first_treated_period_ids(x)]
end

# multiple dispatch for length of pre-treatment periods

function length_T₀(x::BalancedPanelCov{SingleUnitTreatment{Continuous}})
    first_treated_period_ids(x) - 1
end

function length_T₀(x::BalancedPanelCov{MultiUnitTreatment{Simultaneous{Continuous}}})
    first_treated_period_ids(x) .- 1
end

function length_T₀(x::BalancedPanelQ{SingleUnitTreatment{Continuous}})
    first_treated_period_ids(x) - 1
end

function length_T₀(x::BalancedPanelQ{MultiUnitTreatment{Simultaneous{Continuous}}})
    first_treated_period_ids(x) .- 1
end

# multiple dispatch for length of post-treatment periods
function length_T₁(x::BalancedPanelCov{SingleUnitTreatment{Continuous}})
    size(x.Y, 2) .- first_treated_period_ids(x) + 1
end

function length_T₁(x::BalancedPanelCov{MultiUnitTreatment{Simultaneous{Continuous}}})
    size(x.Y[treated_ids(x), :], 2) .- first_treated_period_ids(x) .+ 1
end

function length_T₁(x::BalancedPanelQ{SingleUnitTreatment{Continuous}})
    size(x.Y, 2) .- first_treated_period_ids(x) + 1
end

function length_T₁(x::BalancedPanelQ{MultiUnitTreatment{Simultaneous{Continuous}}})
    size(x.Y[treated_ids(x), :], 2) .- first_treated_period_ids(x) .+ 1
end

# multiple dispatch for get_y₁₀ (the pre-treatment outcomes for the treated units)


#function get_y₁₀(x::BalancedPanelCov{SingleUnitTreatment{Continuous}})
#    x.Y[treated_ids(x), 1:first_treated_period_ids(x)-1]
#end

#function get_y₁₀(x::BalancedPanelQ{SingleUnitTreatment{Continuous}})
#    x.Y[treated_ids(x), 1:first_treated_period_ids(x)-1]
#end


#### big constructors 

## Balanced Panel with covariates stored in a 3rd dimension slice of the outcome array Y 

## DP: Version with covariates 
# Constructor for single continuous treatment - returns BalancedPanel{SingleUnitTreatment{Continuous}}
function BalancedPanelCov(df::DataFrame, 
                        treatment_assignment::Pair{T1, T2};
                        id_var = nothing, 
                        t_var = nothing, 
                        outcome_var = nothing, 
                        covariates::Union{Nothing, Vector{Symbol}, Symbol} = nothing,
                        sort_inplace = false) where T1 where T2 <: Union{Date, Int}
    
    if isnothing(covariates)    # fallback to BalancedPanel constructor if no covariates are specified
        return BalancedPanel(df, treatment_assignment; id_var = id_var, t_var = t_var, outcome_var = outcome_var, sort_inplace = sort_inplace)
    elseif typeof(covariates) == Symbol
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

    # Outcome matrix with covariates
    Y = zeros(panel_y_eltype(df, outcome_var, covariates), (size(W)..., length(covariates)+1))
    for (row, i) ∈ enumerate(is), (col, t) ∈ enumerate(ts)
        Y[row, col,1] = only(df[(df[!, id_var] .== i) .& (df[!, t_var] .== t), outcome_var])
    end
    # add covariates in 3rd dimension
    for covariate in eachindex(covariates)
        for (row, i) ∈ enumerate(is), (col, t) ∈ enumerate(ts)
            Y[row, col,1+covariate] = only(df[(df[!, id_var] .== i) .& (df[!, t_var] .== t), covariates[covariate]])
        end
    end

    return BalancedPanelCov{SingleUnitTreatment{Continuous}}(W, Y, df, id_var, t_var, outcome_var, ts, is)  
end

### multi unit treatment with covariates 
# DP Fallback method with covariates - if the length of treatment assignment is one use single treatment method above
function BalancedPanelCov(df::DataFrame, treatment_assignment::AbstractVector{<:Pair};
                                id_var = nothing, 
                                t_var = nothing, 
                                outcome_var = nothing, 
                                covariates::Union{Nothing, Vector{Symbol}, Symbol} = nothing,
                                sort_inplace = false)  

    # fallback singlue unit treatment
    if (typeof(treatment_assignment) <: Vector && length(treatment_assignment) == 1)
        return BalancedPanelCov(df, only(treatment_assignment); id_var = id_var, t_var = t_var,
                                outcome_var = outcome_var, covariates = covariates, sort_inplace = sort_inplace)
    end

    # fallback no covariates 
    if isnothing(covariates)
        return BalancedPanel(df, treatment_assignment; id_var = id_var, t_var = t_var, outcome_var = outcome_var, sort_inplace = sort_inplace)
    elseif typeof(covariates) == Symbol
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
    Y = zeros(panel_y_eltype(df, outcome_var, covariates), (size(W)..., length(covariates)+1))
    for (row, i) ∈ enumerate(is), (col, t) ∈ enumerate(ts)
        Y[row, col,1] = only(df[(df[!, id_var] .== i) .& (df[!, t_var] .== t), outcome_var])
    end
    # add covariates in 3rd dimension
    for covariate in eachindex(covariates)
        for (row, i) ∈ enumerate(is), (col, t) ∈ enumerate(ts)
            Y[row, col,1+covariate] = only(df[(df[!, id_var] .== i) .& (df[!, t_var] .== t), covariates[covariate]])
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

    BalancedPanelCov{MultiUnitTreatment{uttype{tdtype}}}(W, Y, df, id_var, t_var, outcome_var, ts, is)  
end

### now also add instances with baseline weights Q (either a dyadic tensor object or an N × T matrix) 

"""
    Function to build the 3-dimensional tensor of baseline weights.
    Note that currently doesn't work in pipeline because the 
    cohort selection leads to df being a subset of i,j and t of those in df_q !
    (need to implement a filter first)
"""
function assemble_tensor!(Q, df_q, is, ts, id_var, q_id_var, t_var, q_var)#; combine=:overwrite)
    # Build compact ID→index maps (small Dicts, used many times)
    i_idx = Dict{eltype(is), Int}(zip(is, eachindex(is)))
    j_idx = Dict{eltype(is), Int}(zip(is, eachindex(is)))  # if same id space; otherwise use q_ids
    t_idx = Dict{eltype(ts), Int}(zip(ts, eachindex(ts)))

    # Local column refs (faster than row[:col] / property access in a tight loop)
    id_col  = df_q[!, id_var]
    qid_col = df_q[!, q_id_var]
    t_col   = df_q[!, t_var]
    v_col   = df_q[!, q_var]

    @inbounds for r in eachindex(id_col)
        i = i_idx[id_col[r]]       # KeyError if unseen; handle if needed
        j = j_idx[qid_col[r]]
        t = t_idx[t_col[r]]

        #if combine === :overwrite
        Q[i, j, t] = v_col[r]
        #elseif combine === :sum
        #    Q[i, j, t] += v_col[r]
        #elseif combine === :max
        #    Q[i, j, t] = max(Q[i, j, t], v_col[r])
        #else
        #    Q[i, j, t] = v_col[r]  # default
        #end
    end
    return Q
end

### PanelMakers for BalancedPanelQ
# Constructor for single continuous treatment - returns BalancedPanel{SingleUnitTreatment{Continuous}}
function BalancedPanelQ(df::DataFrame, 
                        treatment_assignment::Pair{<:Any, <:Union{Date, Int}},
                        df_q::Union{Nothing, DataFrame} = nothing;
                        id_var = nothing, 
                        t_var = nothing, 
                        outcome_var = nothing, 
                        q_var = nothing,                # the variable name of the weights in either df or df_q
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
        Y = zeros(panel_y_eltype(df, outcome_var, covariates), (size(W)..., length(covariates)+1))
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

        #assemble_tensor!(Q, df_q, is, ts, id_var, q_id_var, t_var, q_var)
        if string(t_var) ∈ string.(names(df_q))      # this controls whether weights in df_q are time-variant
            q_lookup = Dict{Tuple{eltype(is), eltype(is), eltype(ts)}, eltype(df_q[!, q_var])}()

            for row in eachrow(df_q)
                key = (row[id_var], row[q_id_var], row[t_var])
                q_lookup[key] = row[q_var]  # Overwrites duplicates, change if needed
            end

            # Now loop and fill Q efficiently
            for (row, i) ∈ enumerate(is), (col, t) ∈ enumerate(ts), (slice, j) ∈ enumerate(is)
                Q[row, col, slice] = get(q_lookup, (i, j, t), 0)
            end
        else
            q_lookup = Dict{Tuple{eltype(is), eltype(is)}, eltype(df_q[!, q_var])}()

            for row in eachrow(df_q)
                key = (row[id_var], row[q_id_var])
                q_lookup[key] = row[q_var]  # Overwrites duplicates, change if needed
            end

            # Now loop and fill Q efficiently, note here we only look up i and j and fill every t with this
            for (row, i) ∈ enumerate(is), (col, t) ∈ enumerate(ts), (slice, j) ∈ enumerate(is)
                Q[row, col, slice] = get(q_lookup, (i, j), 0)
            end
        end
        # use of get() ensures dyads not present in df get assigned a baseline weight of 0

        #for (row, i) ∈ enumerate(is), (col, t) ∈ enumerate(ts), (slice, j) ∈ enumerate(is)
        #    Q[row, col, slice] = only(df_q[(df_q[!, id_var] .== i) .& (df_q[!, q_id_var] .== j) .& (df_q[!, t_var] .== t), q_var])
        #end
        
    else
        if isnothing(df_q) && !isnothing(q_var)
            df_q = select(df, [id_var, t_var, q_var])
        elseif isnothing(df_q) && isnothing(q_var)
            error(ArgumentError("Neither dataframe for weights df_q provided nor variable q_var for weights specified."))
        end

        if string(t_var) ∉ string.(names(df_q)) && length(unique(df_q[!,id_var])) == length(df_q[!,id_var])
            # case where weights are static
            Q = zeros(eltype(df_q[!, q_var]), (size(W, 1)))
            q_lookup = Dict{eltype(is), eltype(df_q[!, q_var])}()

            for row in eachrow(df_q)
                q_lookup[row[id_var]] = row[q_var]
            end

            for (row, i) ∈ enumerate(is)
                Q[row] = get(q_lookup, i, zero(eltype(df_q[!, q_var])))
            end
            # cases where q_var is missing for potential control units get assigned a baseline weight 0
        else
            # case where weights change in time
            @assert string(t_var) ∈ string.(names(df_q))
            Q = zeros(eltype(df_q[!, q_var]), (size(W)))
            q_lookup = Dict{Tuple{eltype(is), eltype(ts)}, eltype(df_q[!, q_var])}()

            for row in eachrow(df_q)
                q_lookup[(row[id_var], row[t_var])] = row[q_var]
            end

            for (row, i) ∈ enumerate(is), (col, t) ∈ enumerate(ts)
                Q[row, col] = get(q_lookup, (i, t), zero(eltype(df_q[!, q_var])))
            end
            # cases where q_var is missing for potential control units get assigned a baseline weight 0
        end
    end
    # XXX Problem: Kullback-Leibler Divergence not defined for instances with zero baseline weights. Need to add step to purge controls with zero baselines

    return BalancedPanelQ{SingleUnitTreatment{Continuous}}(W, Y, Q, df, id_var, t_var, outcome_var, ts, is)  
end




# Fallback method - if the length of treatment assignment is one use single treatment method above
function BalancedPanelQ(df::DataFrame, treatment_assignment::AbstractVector{<:Pair},
                                df_q::Union{Nothing, DataFrame} = nothing;
                                id_var = nothing, 
                                t_var = nothing, 
                                outcome_var = nothing, 
                                q_var = nothing,
                                q_id_var = nothing,
                                q_dyadic::Bool = true,
                                covariates::Union{Nothing, Vector{Symbol}, Symbol} = nothing,
                                sort_inplace = false)

    if (typeof(treatment_assignment) <: Vector && length(treatment_assignment) == 1)
        return BalancedPanelQ(df, only(treatment_assignment), df_q; id_var = id_var, t_var = t_var,
                                outcome_var = outcome_var, q_var = q_var, q_id_var = q_id_var, q_dyadic = q_dyadic, covariates = covariates, sort_inplace = sort_inplace)
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
        Y = zeros(panel_y_eltype(df, outcome_var, covariates), (size(W)..., length(covariates)+1))
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

        #assemble_tensor!(Q, df_q, is, ts, id_var, q_id_var, t_var, q_var)
        if string(t_var) ∈ string.(names(df_q))      # this controls whether weights in df_q are time-variant
            q_lookup = Dict{Tuple{eltype(is), eltype(is), eltype(ts)}, eltype(df_q[!, q_var])}()

            for row in eachrow(df_q)
                key = (row[id_var], row[q_id_var], row[t_var])
                q_lookup[key] = row[q_var]  # Overwrites duplicates, change if needed
            end

            # Now loop and fill Q efficiently
            for (row, i) ∈ enumerate(is), (col, t) ∈ enumerate(ts), (slice, j) ∈ enumerate(is)
                Q[row, col, slice] = get(q_lookup, (i, j, t), 0)
            end
        else
            q_lookup = Dict{Tuple{eltype(is), eltype(is)}, eltype(df_q[!, q_var])}()

            for row in eachrow(df_q)
                key = (row[id_var], row[q_id_var])
                q_lookup[key] = row[q_var]  # Overwrites duplicates, change if needed
            end

            # Now loop and fill Q efficiently, note here we only look up i and j and fill every t with this
            for (row, i) ∈ enumerate(is), (col, t) ∈ enumerate(ts), (slice, j) ∈ enumerate(is)
                Q[row, col, slice] = get(q_lookup, (i, j), 0)
            end
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

        if string(t_var) ∉ string.(names(df_q)) && length(unique(df_q[!,id_var])) == length(df_q[!,id_var])
            Q = zeros(eltype(df_q[!, q_var]), (size(W, 1)))
            q_lookup = Dict{eltype(is), eltype(df_q[!, q_var])}()

            for row in eachrow(df_q)
                q_lookup[row[id_var]] = row[q_var]
            end

            for (row, i) ∈ enumerate(is)
                Q[row] = get(q_lookup, i, zero(eltype(df_q[!, q_var])))
            end
        else
            @assert string(t_var) ∈ string.(names(df_q))
            Q = zeros(eltype(df_q[!, q_var]), (size(W)))
            q_lookup = Dict{Tuple{eltype(is), eltype(ts)}, eltype(df_q[!, q_var])}()

            for row in eachrow(df_q)
                q_lookup[(row[id_var], row[t_var])] = row[q_var]
            end

            for (row, i) ∈ enumerate(is), (col, t) ∈ enumerate(ts)
                Q[row, col] = get(q_lookup, (i, t), zero(eltype(df_q[!, q_var])))
            end
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
