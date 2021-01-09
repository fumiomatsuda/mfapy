#-------------------------------------------------------------------------------
# Name:        mfapy Example_1_toymodel.py Example code of mfapy
#-------------------------------------------------------------------------------
import mfapy
# Model construction
reactions, reversible, metabolites, fragments\
 = mfapy.mfapyio.load_metabolic_model("example_1_toymodel_model.txt")
model = mfapy.metabolicmodel.MetabolicModel(reactions, reversible,\
 metabolites, fragments)
# Addition of constraints
state = model.load_states("Example_1_toymodel_status.csv", format = 'csv')
model.set_constraints_from_state_dict(state)
model.update()
# Generation of CarbonSource instance
cs1 = model.generate_carbon_source_templete()
cs1.set_each_isotopomer('AcCoA', {'#10':0.5})
cs2 = model.generate_carbon_source_templete()
cs2.set_each_isotopomer('AcCoA', {'#11':0.5})
# Load MDV data
mdv1 = model.load_mdv_data("Example_1_MDV1.txt")
mdv2 = model.load_mdv_data("Example_1_MDV2.txt")
# Flux estimation Step 1: Setting experimments
model.set_experiment('ex1', mdv1, cs1)
model.set_experiment('ex2', mdv2, cs2)
# Flux estimation step 2: Generation of intical flux vectors
endstate, flux = model.generate_initial_states(50, 4, method ="parallel")
# Flux estimation step 3: Fitting model
for method in ["GN_CRS2_LM", "LN_PRAXIS", "SLSQP"]:
    endstate, RSS, flux = model.fitting_flux(method = method, flux = flux)
model.show_results([("final", flux[0])]) # Show result
# Estimation of 95% CI
ci_edge = model.generate_ci_templete(targets = [('reaction', "v3")])
ci = model.search_ci(ci_edge, flux[0], method = 'grid')
lb = ci['data'][('reaction', "v3")]['lower_boundary']
ub = ci['data'][('reaction', "v3")]['upper_boundary']
print("v3", "Lower bondary:",lb, "Upper boundary:", ub)
