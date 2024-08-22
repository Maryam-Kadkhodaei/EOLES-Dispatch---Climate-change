"""
EOLES Dispatch 30 year test_add_climate
"""

"""
Eoles Model from Behrang Shirizadeh, Quentin Perrier and Philippe Quirion, May 2021
Written in Python by Nilam De Oliveira-Gill, June 2021
"""

"""IMPORTS

Import modules and libraries needed for the programm 
"""
import pyomo.environ as pyo
from pyomo.opt import SolverFactory
import pandas as pd
import csv
import time
import sys

#Initialize time to measure the execution time
start_time = time.time()

#We set the name of the model here, it will be used in outputs name
model_name = ""

"""INITIALISATION OF THE MODEL"""

model = pyo.ConcreteModel()

#Dual Variable, used to get the marginal value of an equation.
model.dual = pyo.Suffix(direction=pyo.Suffix.IMPORT)

"""INPUTS"""
time_step=2
num_of_years =30

    # Installed capacity
#capa = pd.read_csv(r'D:\input for dispatch run\gen_capa_CNRM-CERFACS-CNRM-CM5_CNRM_historical.csv', index_col=[0], header=None)
#capa = capa.squeeze().copy()

    # volume of storage
#capacity = pd.read_csv(r'D:\input for dispatch run\str_charge_energy_CNRM-CERFACS-CNRM-CM5_CNRM_historical.csv', index_col=[0], header=None)
#capacity = capacity.squeeze().copy()

    # storing capacity
#s = pd.read_csv(r'D:\input for dispatch run\str_charge_power_CNRM-CERFACS-CNRM-CM5_CNRM_historical.csv', index_col=[0], header=None)
#s = s.squeeze().copy()

    # Production profiles of VRE
#load_factor = pd.read_csv("inputs_30y/vre_CNRM-CERFACS-CNRM-CM5_CNRM_historical_resampled_2h.csv", index_col=[0, 1], header=None)
#load_factor = load_factor.squeeze().copy()

    # Demand profile in each our in GW
#demand = pd.read_csv("inputs_30y/demand_30y_CNRM-CERFACS-CNRM-CM5_CNRM_historical_calibrated_2h_ts.csv",index_col=0,  header=None)
#demand = demand.squeeze().copy()

    # Monthly lake inflows in GWh
lake_inflows = pd.read_csv("inputs_30y/lake_inflow_30years.csv", index_col=0, header=None)
lake_inflows = lake_inflows.squeeze().copy()


    # Existing capacities of the technologies by December 2017 in GW
capa_ex = pd.read_csv("inputs_30y/existing_capas_elec_new.csv", index_col=0, header=None)
capa_ex = capa_ex.squeeze().copy()

    # Maximum capacities of the technologies in GW
capa_max = pd.read_csv("inputs_30y/max_capas_elec_new.csv", index_col=0, header=None)
capa_max = capa_max.squeeze().copy()

    # Fixed capacities of the technologies in GW
fix_capa = pd.read_csv("inputs_30y/fix_capas.csv", index_col=0,  header=None)
fix_capa = fix_capa.squeeze().copy()

    # Annualized power capex cost in M€/GW/year
capex = pd.read_csv("inputs_30y/annuities_elec_new.csv", index_col=0,  header=None)
capex = capex.squeeze().copy()

    # Annualized energy capex cost of storage technologies in M€/GWh/year
capex_en = pd.read_csv("inputs_30y/str_annuities_elec_new.csv", index_col=0, header=None)
capex_en = capex_en.squeeze().copy()

    # Annualized fixed operation and maintenance costs M€/GW/year
fOM = pd.read_csv("inputs_30y/fO&M_elec_new.csv", index_col=0, header=None)
fOM = fOM.squeeze().copy()

    # Variable operation and maintenance costs in M€/GWh
vOM = pd.read_csv("inputs_30y/vO&M_elec_new.csv", index_col=0,  header=None)
vOM = vOM.squeeze().copy()

    # Charging related annuity of storage in M€/GW/year
s_capex = pd.read_csv("inputs_30y/s_capex.csv", index_col=0,  header=None)
s_capex = s_capex.squeeze().copy()

    # Charging related fOM of storage in M€/GW/year
s_opex = pd.read_csv("inputs_30y/s_opex.csv", index_col=0,  header=None)
s_opex = s_opex.squeeze().copy()

    # Charging efficiency of storage technologies
eta_in = pd.read_csv("inputs_30y/eta_in.csv", index_col=0,  header=None)
eta_in = eta_in.squeeze().copy()

    # Discharging efficiency of storage technolgoies
eta_out = pd.read_csv("inputs_30y/eta_out.csv", index_col=0, header=None)
eta_out = eta_out.squeeze().copy()

    # Existing storage capacity in GWh
capacity_ex = pd.read_csv("inputs_30y/capacity_ex.csv", index_col=0,  header=None)
capacity_ex = capacity_ex.squeeze().copy()
    #Parameters in miscellaneous.csv :
        # eta_ccgt                  : efficiency of CCGT power plants with CCS
        # delta                     : load variation factor
        # H2_demand                 : hourly hydrogen demand on top of the storage
        # eta_electrolysis          : efficiency of Electrolysis
        # phs_discharging_lower     : lower bounds for capa(phs)
        # phs_discharging_upper     : upper bounds for capa(phs)
        # phs_charging_lower        : lower bounds for s(phs) 
        # phs_charging_upper        : upper bounds for s(phs)
        # phs_energy_lower          : lower bounds for capacity(phs)
        # phs_energy_upper          : upper bounds for capacity(phs)
        # first_month               : first month of demand
        # first_year                : first year of demand

miscellaneous = pd.read_csv("inputs_30y/miscellaneous.csv",index_col=0, header=None)
miscellaneous = miscellaneous.squeeze().copy()
"""SET HOUR BY MONTHS

Take the number of hour in the demand file.
Set up hours per months."""

capa_dict = capa.to_dict().copy()
capacity_dict = capacity.to_dict().copy()
charge_power_capa_dict = s.to_dict().copy()

#Maximum generation bounds
def generation_bounds(model,i,h):
        return (None, capa_dict[i]*time_step)

#Maximum stored energy bound
def str_soc_bounds(model,i,h):
        return (None, capacity_dict[i])

#charging power capacity of the battery
def str_ch_p_bounds(model,i,h):
        return (None, charge_power_capa_dict[i]*time_step)



first_year = miscellaneous['first_year']
first_month = miscellaneous['first_month']
first_hour = 0
last_hour = len(demand)
days_in_feb = 672


hours_by_months = {1: 744/time_step, 2: days_in_feb/time_step, 3: 744/time_step, 4: 720/time_step, 5: 744/time_step, 6: 720/time_step, 7: 744/time_step, 8: 744/time_step, 9: 720/time_step, 10: 744/time_step, 11: 720/time_step,
                   12: 744/time_step}

months_count = 12*num_of_years
i = 1
j = first_month
hour = 0
months_hours={}
while i <= months_count:
    months_hours[i] = range(int(hour), int(hour + hours_by_months[j]))
    hour += hours_by_months[j]
    j += 1
    i += 1
    if j == 13:
        j = 1

"""SETS

Definition of set as an object of the model
"""

#Range of hour in one year
model.h = \
    pyo.RangeSet(first_hour,last_hour-1)
#Months
model.months = \
    pyo.RangeSet(1,months_count)
#Technologies
model.tec = \
    pyo.Set(initialize=["onshore", "pv_g", "pv_c", "river", "lake",  "biogas2",  "ccgt", "nuc", "h2_ccgt", "phs",  "battery4", "electrolysis", "hydrogen"])
#Power plants
model.gen = \
    pyo.Set(initialize=[ "onshore", "pv_g", "pv_c", "river", "lake","ccgt","nuc"])
#Variables Technologies
model.vre = \
    pyo.Set(initialize=[ "onshore", "pv_g", "pv_c", "river"])
#
model.balance = \
    pyo.Set(initialize=[ "lake","nuc","phs","battery4","h2_ccgt","ccgt"])
#Storage Technologies
model.str = \
    pyo.Set(initialize=["phs", "battery4", "hydrogen"])
#Storage Technologies
model.str_noH2 = \
    pyo.Set(initialize=["phs", "battery4"])
#Battery Storage
model.battery = \
    pyo.Set(initialize=["battery4"])

"""PARAMETERS"""

#Set the hydrogen demand for each hour
H2_demand = {}
for hour in model.h:
    H2_demand[hour] = miscellaneous['H2_demand']


"""BOUNDS VALUES

Set initial value for variables.
There is a function for each variable with bounds.
The function return the lower and the upper value.
"""


"""VARIABLES

Definition of variable as an object of the model
"""

    # Hourly energy generation in GWh/h
model.gene = \
    pyo.Var(((tec, h) for tec in model.tec for h in model.h), within=pyo.NonNegativeReals,initialize=0, bounds=generation_bounds)


    # Hourly lost load generation in GWh/h
model.ll = \
    pyo.Var(( h for h in model.h), within=pyo.NonNegativeReals,initialize=0)

    # Hourly electricity input of battery storage GW
model.storage = \
    pyo.Var(((storage, h) for storage in model.str for h in model.h), within=pyo.NonNegativeReals,initialize=0,bounds=str_ch_p_bounds)

    # Energy stored in each storage technology in GWh = Stage of charge
model.stored = \
    pyo.Var(((storage, h) for storage in model.str for h in model.h), within=pyo.NonNegativeReals,initialize=0, bounds=str_soc_bounds)





"""CONSTRAINTS RULE

Set up a function which will return the equation of the constraint.
"""



def combustion_2_constraint_rule(model, h):
    """Get constraint on the relationship of combustible technologies"""

    return model.gene['ccgt', h] == model.gene['biogas2',h]*miscellaneous['eta_ccgt']


def storing_constraint_rule(model, h, storage_tecs):
    """Get constraint on storing."""

    hPOne = h+1 if h<(last_hour-1) else 0
    charge = model.storage[storage_tecs, h] * eta_in[storage_tecs]
    discharge =  model.gene[storage_tecs, h] / eta_out[storage_tecs]
    flux = charge - discharge
    return model.stored[storage_tecs, hPOne] == model.stored[storage_tecs, h] + flux

def storage_constraint_rule(model,storage_tecs):
    """Get constraint on stored energy to be equal at the end than at the start."""

    first = model.stored[storage_tecs, first_hour]
    last = model.stored[storage_tecs, last_hour-1]
    charge = model.storage[storage_tecs, last_hour-1] * eta_in[storage_tecs]
    discharge = model.gene[storage_tecs, last_hour-1] / eta_out[storage_tecs]
    flux = charge - discharge
    return first == last + flux

def lake_reserve_constraint_rule(model, month):
    """Get constraint on maximum monthly lake generation."""

    return sum(model.gene['lake', hour] for hour in months_hours[month]) <= lake_inflows[month] * 1000


def biogas_constraint_rule(model):
    """Get constraint on biogas."""

    gene_biogas = sum(model.gene['biogas2',hour] for hour in model.h)

    return gene_biogas <= miscellaneous['max_biogas'] * 1000*num_of_years

def hydrogen_balance_constraint_rule(model,h):
    """Get constraint on hydrogen's balance."""

    gene_e_h = model.gene['electrolysis',h]+model.gene['hydrogen',h]
    dem_sto = model.gene['h2_ccgt',h]/miscellaneous['eta_h2_ccgt'] + H2_demand[h] + model.storage['hydrogen',h]
    return gene_e_h == dem_sto


def adequacy_constraint_rule(model, h):
    """Get constraint for 'supply/demand relation'"""
    vre = time_step*sum(capa[vre] * load_factor[vre, h] for vre in model.vre) 
    sto = sum(model.storage[str_noH2, h] for str_noH2 in model.str_noH2)
    gene_electrolysis = model.gene['electrolysis',h] / miscellaneous['eta_electrolysis']
    return vre + sum(model.gene[balance, h] for balance in model.balance)+model.ll[h] >= (demand[h] + sto + gene_electrolysis)

def objective_rule(model):
    """Get constraint for the final objective function."""
    fixed_cost= (sum((capa[tec] - capa_ex[tec]) * capex[tec] for tec in model.tec) \
           + sum((capacity[storage_tecs]-capacity_ex[storage_tecs]) * capex_en[storage_tecs] for storage_tecs in model.str)\
           + sum(capa[tec] * fOM[tec] for tec in model.tec)\
           + sum(s[storage_tecs] * (s_opex[storage_tecs] + s_capex[storage_tecs]) for storage_tecs in model.str))/1000
    
    variable_cost= sum(sum(model.gene[tec, h] * vOM[tec] for h in model.h) for tec in model.tec)/1000
    ll_cost = sum(model.ll[h] * 10 for h in model.h)/1000

    return fixed_cost*num_of_years + variable_cost + ll_cost

"""CONSTRAINT CREATION

Create the constraint as an object of the model with the function declared earlier as a rule.
"""

model.combustion_2_constraint = \
    pyo.Constraint(model.h, rule=combustion_2_constraint_rule)

model.storing_constraint = \
    pyo.Constraint(model.h,model.str, rule=storing_constraint_rule)

model.storage_constraint = \
    pyo.Constraint(model.str, rule=storage_constraint_rule)

model.lake_reserve_constraint = \
    pyo.Constraint(model.months, rule=lake_reserve_constraint_rule)

model.biogas_constraint = \
    pyo.Constraint(rule=biogas_constraint_rule)

model.hydrogen_balance_contraint = \
    pyo.Constraint(model.h,rule=hydrogen_balance_constraint_rule)

model.adequacy_constraint = \
    pyo.Constraint(model.h, rule=adequacy_constraint_rule)

#Creation of the objective -> Cost
model.objective = pyo.Objective(rule=objective_rule)

"""SOLVE STATEMENT

Choice of the solver.
You can remove the '#' in the third line to display the output of the solver.
"""

opt = SolverFactory('gurobi')
opt.options['FeasibilityTol'] = tol
opt.options['OptimalityTol'] = tol

results = opt.solve(model)


#model.display()
fixed_cost = ((sum((capa[tec] - capa_ex[tec]) * capex[tec] for tec in model.tec) \
           + sum((capacity[storage_tecs]-capacity_ex[storage_tecs]) * capex_en[storage_tecs] for storage_tecs in model.str)\
           + sum(capa[tec] * fOM[tec] for tec in model.tec)\
           + sum(s[storage_tecs] * (s_opex[storage_tecs] + s_capex[storage_tecs]) for storage_tecs in model.str))/1000)*num_of_years

variable_cost =     sum(sum(model.gene[tec, h] * vOM[tec] for h in model.h) for tec in model.tec)/1000

ll_cost = sum(model.ll[h] * 10 for h in model.h)/1000
"""SET OUTPUTS VARIABLES"""

#Dictionnary which will set a little definition for each technology in the model.
technologies_definition = {
    "onshore" : "onshore wind",
    "pv_g" : "pv grounded",
    "pv_c" : "pv commercial",
    "river" : "run-of-river hydro",
    "lake" : "lake and reservoirs",
    "biogas2" : "biogas for ccgt",
    "ccgt" : "combined cycle gas turbine",
    "nuc" : "nuclear",
    "h2_ccgt" : "combined cycle gas turbine using hydrogen",
    "phs" : "pumped hydroelectric energy storage",
    "battery4" : "4 hours battery",
    "electrolysis" : "electrolysis",
    "hydrogen" : "hydrogen removed from storage",
}

    # The whole demand per year in TWh
sumdemand = sum(demand[hour] for hour in model.h) / 1000
    # The whole electricity demand for hydrogen per year in TWh
dem_hydrogen = sum(H2_demand[hour] for hour in model.h) / 1000
    # The whole generation per year in TWh
sumgene = sum(pyo.value(model.gene[gen,hour]) for hour in model.h for gen in model.gen) / 1000
    # The whole energy budgeted for reserves per year in TWh

    # Overall yearly energy generated by the technology in TWh
gene_tec = {}
for tec in model.tec:
    gene_tec[tec] = sum(pyo.value(model.gene[tec,hour]) for hour in model.h) / 1000

    # The whole electricity input for storage per year in TWh
nSTORAGE = {}
for storage in model.str:
    for hour in model.h:
        nSTORAGE[(storage,hour)] = pyo.value(model.storage[storage,hour])

    # Electricity cost per MWh produced (euros/MWh)
lcoe_sys1 = pyo.value(model.objective)*1000/sumgene

    # Yearly storage related loss in % of power production and in TWh
str_loss_percent = 100*(sum(pyo.value(model.storage[storage,hour]) for storage in model.str for hour in model.h)-\
sum(gene_tec[storage]*1000 for storage in model.str)) / (sumgene*1000)
str_loss_TWh = gene_tec['electrolysis']/miscellaneous['eta_electrolysis'] - dem_hydrogen/miscellaneous['eta_electrolysis'] - gene_tec['h2_ccgt']
for storage in model.str:
    if storage != 'hydrogen' :
        str_loss_TWh += sum(nSTORAGE[storage,hour] for hour in model.h)/1000 - gene_tec[storage]

    # Load curtailment in % of power production and in TWh
lc_percent = (100*(sumgene - sumdemand - dem_hydrogen/0.75)/sumgene) - str_loss_percent
lc_TWh = (sumgene - sumdemand - dem_hydrogen/0.75) - str_loss_TWh

    # Dual values
spot_price = {}
gas_price2 = {}
for hour in model.h:
    spot_price[hour] = - 1000000 * model.dual[model.adequacy_constraint[hour]]
    gas_price2[hour] = -1000000 * model.dual[model.combustion_2_constraint[hour]]




    # Marginal Cost
marginal_cost = sum(spot_price[hour] for hour in model.h) / (last_hour)

    # Average cost of hydrogen (euros/kg)
lcoh_1 = pyo.value(capa['electrolysis'])*(capex['electrolysis']+fOM['electrolysis'])
lcoh_2 =  sum(pyo.value(model.gene['electrolysis',hour])*(vOM['electrolysis']+\
    (spot_price[hour]/1000)) for hour in model.h)
lcoh_3 = capex_en['hydrogen']*pyo.value(capacity['hydrogen'])
lcoh_4 = sum(pyo.value(model.gene['electrolysis',hour]) for hour in model.h)
lcoh = (lcoh_1 + lcoh_2 + lcoh_3)*33.33 / lcoh_4
    # Electricity cost per MWh consumed (euros/MWh)
lcoe_sys2 = (pyo.value(model.objective)-(lcoh*dem_hydrogen/33.33))*1000/sumdemand


"""OUTPUTS
    There is 4 output files :
        - Summary           : A little summary with the cost and some others data
        - Hourly-Generation : Hourly data
        - Elec_Balance      : Electric Production and Consumption
        - Capacities        : List of capacities by technologies

The try, except loop is here the manage error in the creation of the outputs.
"""

#Summary

cost_dict= {}

cost_dict['total cost bEuro'] = pyo.value(model.objective) / num_of_years
cost_dict['variable cost bEuro'] = pyo.value(variable_cost)/num_of_years
cost_dict['fixed cost bEuro'] = pyo.value(fixed_cost)/num_of_years
cost_dict['ll cost'] = pyo.value(ll_cost)/num_of_years
cost_dict['lcoh'] = lcoh
cost_dict['lcoe_sys1'] = lcoe_sys1

dct1 = {k:[v] for k,v in cost_dict.items()}  # WORKAROUND
df_cost= pd.DataFrame(dct1)
df_cost.to_csv(r'ouputs_30y _correct_vre_gen/cost_{}_{}.csv'.format(model_name,scenario))

#Hourly_Generation
df_generation = pd.DataFrame()
for tec in model.tec:
 df_generation[tec] = pyo.value(model.gene[tec,:])

df_generation['ll'] = pyo.value(model.ll[:])

df_generation['onshore'] = capa['onshore'] * load_factor['onshore', :]*time_step
df_generation['pv_g'] = capa['pv_g'] * load_factor['pv_g', :]*time_step
df_generation['pv_c'] =  capa['pv_c'] * load_factor['pv_c', :]*time_step
df_generation['river'] = capa['river'] * load_factor['river', :]*time_step
df_generation['demand'] =  demand


for str in model.str:
  df_generation[str+'soc'] = pyo.value(model.stored[str, :])

for str in model.str:
  df_generation[str+'ch p '] = pyo.value(model.storage[str, :])

df_generation.to_csv(r"ouputs_30y/generation_{}_{}.csv".format(model_name, scenario))


# price
df_price = pd.DataFrame(spot_price.values())
df_price.to_csv(r"ouputs_30y/price_{}_{}.csv".format(model_name, scenario))

#Elec_balance
balance_dict= {}

for tec in model.tec:
    balance_dict[tec] = sum (pyo.value(model.gene[tec,:]))/1000


dct1 = {k:[v] for k,v in balance_dict.items()}
df_balance= pd.DataFrame(dct1)
df_balance.to_csv(r'ouputs_30y/balance_{}_{}.csv'.format(model_name, scenario))

mult_dict= {}

mult_dict['biogas mult'] = -1000000*model.dual[biogas_constraint_rule[:]]
print(-1000000*model.dual[biogas_constraint_rule])
dct1 = {k:[v] for k,v in mult_dict.items()}
df_mult= pd.DataFrame(dct1)
df_mult.to_csv(r'ouputs_30y _correct_vre_gen/mult_{}_{}.csv'.format(model_name,scenario))
