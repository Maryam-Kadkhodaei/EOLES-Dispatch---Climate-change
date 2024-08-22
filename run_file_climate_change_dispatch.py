# -*- coding: utf-8 -*-
"""
Created on Sun May  7 03:17:59 2023

@author: Maryam
"""
import pandas as pd


#n = 'Climate change' 
#m = '30 years'

model_name_list = ['ICHEC-EC-EARTH_ETH','CNRM-CERFACS-CNRM-CM5_CNRM','CNRM-CERFACS-CNRM-CM5_ETH','MOHC-HadGEM2-ES_CNRM','MPI-M-MPI-ESM-LR_ETH']


senario_dict = {'historical':{'demand':'historical','cf':'historical'},
                'rcp85':{'demand':'rcp85','cf':'rcp85'},
                'd_r_cf_h':{'demand':'rcp85','cf':'historical'},
                'd_h_cf_r':{'demand':'historical','cf':'rcp85'}}

scenario_list = ['historical']
tol = 1e-6


sim_type = 'D climate change 30_y'

for scenario in  scenario_list:
      for model_name in model_name_list:


          # Installed capacity
          capa = pd.read_csv(r'inputs_30y/capa_{}_{}.csv'.format(model_name, scenario), index_col=[0], header=None)
          capa = capa.squeeze().copy()
          print(capa)

          # volume of storage

          capacity =  pd.read_csv(r'inputs_30y/capacity_{}_{}.csv'.format(model_name, scenario), index_col=[0], header=None)
          capacity = capacity.squeeze().copy()
          print(capacity)

          # storing capacity
          s = pd.read_csv(r'inputs_30y/s_{}_{}.csv'.format(model_name, scenario), index_col=[0], header=None)
          s = s.squeeze().copy()


          demand = pd.read_csv(r"inputs_30y/demand_30y_{}_{}_calibrated_2h_ts.csv".format(model_name, senario_dict[scenario]['demand']), index_col=0 , header=None)
          demand = demand.squeeze().copy()


          load_factor = pd.read_csv(r"inputs_30y/vre_{}_{}_resampled_2h.csv".format(model_name, senario_dict[scenario]['cf']),index_col= [0, 1], header=None)
          load_factor = load_factor.squeeze().copy()

          print(sim_type, model_name, scenario)
          exec(open('Eoles_elec_dispatch_30y_2ts.py').read())



