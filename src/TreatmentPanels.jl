module TreatmentPanels_Daniel

    include("TreatmentPanel.jl")
    include("show_and_plot.jl")

    # Export main types
    export TreatmentPanel
    export BalancedPanel 
    export BalancedPanel_maker
    # export UnbalancedPanel - Not yet implemented

    # Export treatment description types
    export TreatmentType
    export SingleUnitTreatment, MultiUnitTreatment
    export TreatmentTimingType
    export Staggered, Simultaneous
    export TreatmentDurationType 
    export Continuous, Discontinuous

    # Export utility functions
    export treated_ids, treated_labels, first_treated_period_ids, first_treated_period_labels, length_T₀, length_T₁
    export get_y₀₀, get_y₀₁, get_y₁₀, get_y₁₁, decompose_y

    export treatment_periods, construct_W, check_id_t_outcome

end



#=
import TreatmentPanels.TreatmentPanel
import TreatmentPanels.BalancedPanel 
# import TreatmentPanels.UnbalancedPanel - Not yet implemented

# import TreatmentPanels.treatment description types
import TreatmentPanels.TreatmentType
import TreatmentPanels.SingleUnitTreatment, TreatmentPanels.MultiUnitTreatment
import TreatmentPanels.TreatmentTimingType
import TreatmentPanels.Simultaneous
import TreatmentPanels.TreatmentDurationType 
import TreatmentPanels.Continuous

# Export utility functions
import TreatmentPanels.treated_ids, TreatmentPanels.treated_labels, TreatmentPanels.first_treated_period_ids, TreatmentPanels.first_treated_period_labels, TreatmentPanels.length_T₀, TreatmentPanels.length_T₁
import TreatmentPanels.get_y₀₀, TreatmentPanels.get_y₀₁, TreatmentPanels.get_y₁₀, TreatmentPanels.get_y₁₁, TreatmentPanels.decompose_y
=#

#=
import TreatmentPanels_Daniel.TreatmentPanel
import TreatmentPanels_Daniel.BalancedPanel 
import TreatmentPanels_Daniel.BalancedPanel_maker 
# import TreatmentPanels_Daniel.UnbalancedPanel - Not yet implemented

# import TreatmentPanels_Daniel.treatment description types
import TreatmentPanels_Daniel.TreatmentType
import TreatmentPanels_Daniel.SingleUnitTreatment, TreatmentPanels_Daniel.MultiUnitTreatment
import TreatmentPanels_Daniel.TreatmentTimingType
import TreatmentPanels_Daniel.Simultaneous
import TreatmentPanels_Daniel.TreatmentDurationType 
import TreatmentPanels_Daniel.Continuous

# Export utility functions
import TreatmentPanels_Daniel.treated_ids, TreatmentPanels_Daniel.treated_labels, TreatmentPanels_Daniel.first_treated_period_ids, TreatmentPanels_Daniel.first_treated_period_labels, TreatmentPanels_Daniel.length_T₀, TreatmentPanels_Daniel.length_T₁
import TreatmentPanels_Daniel.get_y₀₀, TreatmentPanels_Daniel.get_y₀₁, TreatmentPanels_Daniel.get_y₁₀, TreatmentPanels_Daniel.get_y₁₁, TreatmentPanels_Daniel.decompose_y
=#