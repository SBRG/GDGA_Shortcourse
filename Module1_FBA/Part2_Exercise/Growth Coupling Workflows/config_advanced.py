"""
gcOpt config file
For more optional attributes refer to the default config files in the growth coupling suite directories
"""

import os


# %% mandatory gcOpt configurations
output_dir = os.getcwd() + "/Results/Butanediol_Advanced" # directory in which results are saved
output_file_id = "i8_k7_a7"

biomass_id = "BIOMASS_Ec_iML1515_core_75p37M" # specify the biomass equation ID
growth_rate_fix = 0.1   # fix biomass value


# %% main gcOpt configurations

# solver specific
time_limit = 300 # runtime limit of (gurobi) solver [sec]
processes = 7 # number of threads gurobi is using

# specific number of interventions
num_total_interventions = 8 # number of total interventions

num_deletions = 7 # maximum number of reaction deletions
num_addins = 7 # maximum number of heterologous reaction additions
num_cofeeds = 1 # maximum number of additional metabolite cofeeds, only exchanges which are not source targets are cofeeds
num_carbonsources = 1 # maximum number of carbon sources
num_mediareductions = 0 # not yet implemented

# evaluate GPR relations for each design solution
eval_gpr = True



# %% extensive functionalities
# explicit targets or exclusions

# enforce target space
deletion_targets = []
addin_targets = [] # has no effect
cofeed_targets = []
source_targets = ["EX_arab__L_e", "EX_acgam_e", "EX_succ_e", "EX_gal_e", "EX_asp__L_e",
                  "EX_pro__L_e", "EX_ala__D_e", "EX_tre_e", "EX_man_e", "EX_ser__D_e",
                  "EX_sbt__D_e", "EX_glyc_e", "EX_fuc__L_e", "EX_glcur_e",
                  "EX_glyc3p_e", "EX_xyl__D_e", "EX_lac__L_e", "EX_mnl_e",
                  "EX_g6p_e", "EX_mal__L_e", "EX_mal__D_e", "EX_rib__D_e", "EX_rmn_e",
                  "EX_fru_e", "EX_ac_e", "EX_glc__D_e", "EX_malt_e", "EX_melib_e",
                  "EX_thymd_e", "EX_asn__L_e", "EX_akg_e", "EX_lcts_e", "EX_sucr_e",
                  "EX_uri_e", "EX_gln__L_e", "EX_g1p_e", "EX_f6p_e",
                  "EX_malttr_e", "EX_dad_2_e", "EX_adn_e", "EX_fum_e",
                  "EX_ppa_e", "EX_glyclt_e", "EX_ins_e", "EX_ser__L_e",
                  "EX_thr__L_e", "EX_ala__L_e", "EX_acac_e", "EX_acmana_e", "EX_lyx__L_e",
                  "EX_pyr_e", "EX_galur_e", "EX_all__D_e", "EX_acnam_e",
                  "EX_gam_e", "EX_dha_e" ] # growth was checked on BioCyc (Biolog PM1/2)
mediareduction_targets = []

# explicitly exclude from target space
deletion_exclude_list = ["NADH16pp", "NADH17pp", "NADH18pp", "NADTRHD", "ATPS4rpp",
                         "CYTBO3_4pp", "CYTBDpp", "CYTBD2pp"]
addin_exclude_list = []
cofeed_exclude_list = [] 
source_exclude_list = -1
mediareduction_exclude_list = []


# exclusion lists
subsystem_exclude_list = ['Cell Envelope Biosynthesis',
                          'Exchange',
                          'Inorganic Ion Transport and Metabolism',
                          'Lipopolysaccharide Biosynthesis / Recycling',
                          'Murein Biosynthesis',
                          'Murein Recycling',
                          'Transport, Inner Membrane',
                          'Transport, Outer Membrane',
                          'Transport, Outer Membrane Porin',
                          'tRNA Charging',
                          'TRANSPORT, EXTRACELLULAR',
                          'Transport, extracellular',
                          'Secondary transporters',
                          'Transport, Extracellular'
                          ]

# list of exchanges not being touched if source is an optimization variable
exchanges_not_to_knockout = ['EX_co2_e', 'EX_h2o_e', 'EX_h_e', 'EX_fe2_e', 'EX_fe3_c',
                     'EX_na1_e', 'EX_cl_e', 'EX_k_e', 'EX_o2_e', 'EX_nh4_e']


# disregard exchange reactions as cofeed and carbon source target
exchanges_not_to_add = ["EX_co2_e", "EX_h2_e", "EX_h2s_e", "EX_o2_e", "EX_ch4_e", "EX_o2s_e",
                         "EX_h2o2_e", "EX_h2o_e", "EX_h_e", "EX_n2o_e", "EX_no_e", "EX_cynt_e",
                        "EX_meoh_e", "EX_cyan_e", "EX_mso3_e", "EX_mepn_e", "EX_tcynt_e"
                         ]

# explicitly include exchange reactions as cofeed and carbon source target
exchanges_to_keep = []

# maximum carbon content of cofeed or carbon source targets
num_exchange_carbons_threshold = 8

# %% options for processing deletions targets
# Essentiality of reactions (not considered as deletions targets) is based on the original wild-type model
consider_wildtype_essentiality = True



# %% options for processing carbon source and cofeed choice
maximum_source_non_carbon_flux = 10 # [mmol/gDW/h]
maximum_source_non_carbon_mass_flux = 1.80 # [g/gDW/h]
maximum_source_carbon_flux = 60 # [c-mmol/gDW/h]

# specifically define constraints for nitrogen source or cofeed targets, disabled by default
maximum_source_nitrogen_flux = None # [N-mmol/gDW/h]

user_exclude_list = []




# %% options for processing heterologous reactions (addins)


# assess directions of heterologous reaction by thermodynamic and flux variability analysis
directionality_assessment = True

# if thermodynamic assessment of reactions is done, only consider reactions for which
# Gibbs free energy of reaction is accessable
consider_thermodynamic_valid_reactions_only = True

# currency metabolites do not count as high carbon metabolites
currency_mets = {'h2o', 'co2', 'o2', 'h2o2', 'nh4', 'no2', 'no3', 'no', 'h2s',
                 'so3','so4','h','h2','pi','ppi','coa','accoa','ppcoa','aacoa',
                 'butcoa','succoa','atp','gtp','adp','gdp','amp','gmp','nad',
                 'nadp','nadh','nadph','fad','fadh','na1','ahcys','amet','thf',
                 'mlthf', 'q8h2','q8','mql8','mqn8','2dmmql8','2dmmq8'}


num_carbons_threshold = 35



# %% model processing options
# adapt model bounds and coefficients to avoid numerical issues when solving the MILP
improve_model_numerics = True



