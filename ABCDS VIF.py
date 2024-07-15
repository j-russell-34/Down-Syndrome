# -*- coding: utf-8 -*-
"""
Created on Fri Dec 15 16:33:47 2023

@author: russj13
"""

import statsmodels.formula.api as smf
import pandas as pd
import functions
import numpy as np
from statsmodels.stats.outliers_influence import variance_inflation_factor
import statsmodels.api as sm

# study_cohort = input("Name of study cohort? \n")
# feobv_csv = input("Name of csv with FEOBV data? \n")
# demographics_csv = input("name of csv with demographics? \n")
# cbf_csv = input("name of csv with cholinergic basal forebrain vols? \n")



study_cohort = "ABC-DS"
centiloid_csv = "Centiloid_ABCDS.csv"
demographics_csv = "Participant_Demographics_ABCDS.csv"
cbf_csv = "ABCDS_DnSeg_volumes.csv"
age_csv = "Age_at_Event_ABCDS.csv"
tiv_csv = "ABCDS_SAMSEG.csv"



# import csvs
centiloid_df = pd.read_csv(centiloid_csv)
demo_df = pd.read_csv(demographics_csv)
cbf_df = pd.read_csv(cbf_csv)
age_df = pd.read_csv(age_csv)
tiv_df = pd.read_csv(tiv_csv)

# remove hyphens from headers on sclimbic output
#cbf_df = cbf_df.rename(columns={'Left-Basal-Forebrain': 'Left_ChBF', 'Right-Basal-Forebrain': 'Right_ChBF'})


#Remove values with 0 for ChBF
cbf_df = cbf_df.loc[cbf_df.Left_DnSeg!=0]
cbf_df = cbf_df.loc[cbf_df.Right_DnSeg!=0]

# remove non BL scans and visits

centiloid_df = centiloid_df.loc[centiloid_df.event_code=="bl"]
age_df = age_df.loc[age_df.event_code=="bl"]
demo_df = demo_df.loc[demo_df.event_code=="bl"]
centiloid_df = centiloid_df.reset_index(drop=True)

demo_df = demo_df.reset_index(drop=True)

# get subjects who had an amyloid scan
subject_list_rpts = functions.feobvsubs(centiloid_df, col_id="subject_label")
subject_list = []
for i in subject_list_rpts:
    if i not in subject_list:
        subject_list.append(i)
        
# get subjects who have flair data
flair_subjects_rpts = functions.feobvsubs(tiv_df, col_id="SUBJECT")
flair_subjects = []
for i in flair_subjects_rpts:
    if i not in flair_subjects:
        flair_subjects.append(i)
        
#get subject list who have mri
mri_subjects_rpts = functions.feobvsubs(cbf_df, col_id="SUBJECT")
mri_subjects = []
for i in mri_subjects_rpts:
    if i not in mri_subjects:
        mri_subjects.append(i)
        
#Compare 3 subject lists and make final subject list to use to select subjects from dfs
final_list = set(mri_subjects).intersection(flair_subjects, subject_list)

#Use final list to select data from amyloid, ChBF and TIV dfs
# get just subjects with amyloid scans from BFC df:

centiloidsub_cbf_df = functions.select_feobv_sub(subject_list=final_list, df=cbf_df, col_id="SUBJECT")

# get subjects from amyloid df who have mri data

cbf_centiloid_df = functions.select_feobv_sub(subject_list=final_list, df=centiloid_df, col_id="subject_label")

#get TIV subjects that have all other data
final_tiv_df = functions.select_feobv_sub(subject_list=final_list, df=tiv_df, col_id="SUBJECT")

# sort amyvid and CHBF dataframes to ensure order
cbf_centiloid_df = cbf_centiloid_df.sort_values("subject_label")
centiloidsub_cbf_df = centiloidsub_cbf_df.sort_values("SUBJECT")
final_tiv_df = final_tiv_df.sort_values("SUBJECT")

final_tiv_df = final_tiv_df.reset_index(drop=True)

# edit using newly generated final subject list
# select only subjects that have had centiloid scans and mri for demo and age (in subjects list from demographics CSV
mri_centiloid_demo_df = functions.select_feobv_sub(subject_list=final_list, df=demo_df, col_id="subject_label")
mri_centiloid_age_df = functions.select_feobv_sub(subject_list=final_list, df=age_df, col_id="subject_label")

#sort demographics and age to same order as FEOBV and BCF
mri_centiloid_demo_df = mri_centiloid_demo_df.sort_values("subject_label")
mri_centiloid_age_df = mri_centiloid_age_df.sort_values("subject_label")



# average CH BF vols where multiple scans and select just columns with ChBFvol:
centiloidsub_cbf_df = centiloidsub_cbf_df[["SUBJECT", "Left_DnSeg", "Right_DnSeg"]]
centiloidsub_cbf_df = centiloidsub_cbf_df.groupby(["SUBJECT"]).mean().reset_index()

#divide ChBF by TIV 
centiloidsub_cbf_df['Left_DnSeg'] =  centiloidsub_cbf_df['Left_DnSeg'] / final_tiv_df['samseg_sbtiv']
centiloidsub_cbf_df['Right_DnSeg'] =  centiloidsub_cbf_df['Right_DnSeg'] / final_tiv_df['samseg_sbtiv']

# select ages and add to dataframe for correlation matrix

agecolumn = mri_centiloid_age_df.filter(regex='age_at_visit')
tracer_col = cbf_centiloid_df.filter(regex='tracer')
sex = mri_centiloid_demo_df.filter(regex='de_gender')
df_full = pd.concat([centiloidsub_cbf_df.reset_index(drop=True), agecolumn.reset_index(drop=True)], axis=1)


sex = sex.reset_index(drop=True)



# Not multiple ChBFs with Scilimbic. Can use to average across left ant right if wanted later
# calculate the mean of the bilateral ChBF vols generrating 1 value for Ch1,2,3 and one for Ch4 per participant
# df_full["Ch123_vol"] = df_full.loc[:, ["CH123_L_VOL", "CH123_R_VOL"]].mean(axis=1)
# df_full["Ch4_vol"] = df_full.loc[:, ["CH4_L_VOL", "CH4_R_VOL"]].mean(axis=1)


# select column with centiloid values

centiloid_col = cbf_centiloid_df['centiloid_value']

df_full = pd.concat([df_full.reset_index(drop=True), centiloid_col.reset_index(drop=True)], axis=1)

df_full = pd.concat([df_full.reset_index(drop=True), sex.reset_index(drop=True)], axis=1)

df_full = pd.concat([df_full.reset_index(drop=True), tracer_col.reset_index(drop=True)], axis=1)


#drop NaNs for full_db
cleandffull = df_full.dropna(axis=0)


#center variables
cleandffull["CL_centered"] = cleandffull["centiloid_value"] - np.mean(cleandffull["centiloid_value"])
cleandffull["age_centered"] = cleandffull["age_at_visit"] - np.mean(cleandffull["age_at_visit"])


# iterate through Left and right ChBF vols performing a multiple linear regression for each:
for ChBF in cleandffull.columns[1:3]:
    X = sm.add_constant(cleandffull[["CL_centered", "age_centered", "de_gender"]])
    lm = smf.ols(formula=f"{ChBF} ~ CL_centered + age_centered + de_gender", data=cleandffull).fit()
    print(lm.summary())
    
    vif = pd.DataFrame()
    vif["Variable"] = X.columns
    vif["VIF"] = [variance_inflation_factor(X.values, i) for i in range(X.shape[1])]
    
    print("VIF:")
    print(vif)
    