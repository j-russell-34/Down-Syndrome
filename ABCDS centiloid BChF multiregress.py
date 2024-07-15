import statsmodels.formula.api as smf
import pandas as pd
import functions
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import seaborn as sns
import datetime


# study_cohort = input("Name of study cohort? \n")
# feobv_csv = input("Name of csv with FEOBV data? \n")
# demographics_csv = input("name of csv with demographics? \n")
# cbf_csv = input("name of csv with cholinergic basal forebrain vols? \n")



study_cohort = "ABC-DS"
centiloid_csv = "Centiloid_ABCDS.csv"
demographics_csv = "Participant_Demographics_ABCDS.csv"
cbf_csv = "ABCDS_sclimbic_volumes.csv"
age_csv = "Age_at_Event_ABCDS.csv"
tiv_csv = "ABCDS_SAMSEG.csv"



# import csvs
centiloid_df = pd.read_csv(centiloid_csv)
demo_df = pd.read_csv(demographics_csv)
cbf_df = pd.read_csv(cbf_csv)
age_df = pd.read_csv(age_csv)
tiv_df = pd.read_csv(tiv_csv)

# remove hyphens from headers on sclimbic output
cbf_df = cbf_df.rename(columns={'Left-Basal-Forebrain': 'Left_ChBF', 'Right-Basal-Forebrain': 'Right_ChBF'})


#Remove values with 0 for ChBF
cbf_df = cbf_df.loc[cbf_df.Left_ChBF!=0]
cbf_df = cbf_df.loc[cbf_df.Right_ChBF!=0]

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

#  edit using newly generated final subject list
# select only subjects that have had centiloid scans and mri for demo and age (in subjects list from demographics CSV
mri_centiloid_demo_df = functions.select_feobv_sub(subject_list=final_list, df=demo_df, col_id="subject_label")
mri_centiloid_age_df = functions.select_feobv_sub(subject_list=final_list, df=age_df, col_id="subject_label")

#sort demographics and age to same order as FEOBV and BCF
mri_centiloid_demo_df = mri_centiloid_demo_df.sort_values("subject_label")
mri_centiloid_age_df = mri_centiloid_age_df.sort_values("subject_label")



# average CH BF vols where multiple scans and select just columns with ChBFvol:
centiloidsub_cbf_df = centiloidsub_cbf_df[["SUBJECT", "Left_ChBF", "Right_ChBF"]]
centiloidsub_cbf_df = centiloidsub_cbf_df.groupby(["SUBJECT"]).mean().reset_index()

#divide ChBF by TIV 
centiloidsub_cbf_df['Left_ChBF'] =  centiloidsub_cbf_df['Left_ChBF'] / final_tiv_df['samseg_sbtiv']
centiloidsub_cbf_df['Right_ChBF'] =  centiloidsub_cbf_df['Right_ChBF'] / final_tiv_df['samseg_sbtiv']

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



# iterate through Left and right ChBF vols performing a multiple linear regression for each:
for ChBF in cleandffull.columns[1:3]:
    lm = smf.ols(formula=f"{ChBF} ~ centiloid_value * age_at_visit + de_gender", data=cleandffull).fit()
    print(lm.summary())
    
for ChBF in cleandffull.columns[1:3]:
    lm = smf.ols(formula=f"{ChBF} ~ centiloid_value + de_gender", data=cleandffull).fit()
    print(lm.summary())

for ChBF in cleandffull.columns[1:3]:
    lm = smf.ols(formula=f"{ChBF} ~ centiloid_value * age_at_visit", data=cleandffull).fit()
    print(lm.summary())
    
for ChBF in cleandffull.columns[1:3]:
    lm = smf.ols(formula=f"{ChBF} ~ centiloid_value", data=cleandffull).fit()
    print(lm.summary())


# plot simple linear regression +/- other factors



# for ChBF in ChBFlist:
#     fig, ax = plt.subplots(figsize=(14, 10), facecolor='white')
#     ax.scatter(x=cleandffull['centiloid_value'], y=cleandffull[ChBF], s=75)
#     ax.grid(axis="both")
#     ax.set_facecolor('white')
#     ax.spines['bottom'].set_color('k')
#     ax.spines['left'].set_color('k')
#     ax.spines['top'].set_color('k')
#     ax.spines['right'].set_color('k')
#     plt.xticks(fontsize=28)
#     plt.yticks(fontsize=28)
#     plt.xlabel("Global Amyloid (Centiloids)", fontsize=36)
#     plt.ylabel(f"{ChBF} Volume", fontsize=32)

#     # fit with polyfit
#     coef = np.polyfit(cleandffull['centiloid_value'], cleandffull[ChBF], 1)
#     poly1d_fn = np.poly1d(coef)
#     plt.plot(cleandffull['centiloid_value'], cleandffull[ChBF], 'bo', cleandffull['centiloid_value'], poly1d_fn(cleandffull['centiloid_value']), '-k')

#     correlation = np.corrcoef(cleandffull['centiloid_value'], cleandffull[ChBF])[0, 1]
#     pearson = np.round(stats.pearsonr(cleandffull['centiloid_value'], cleandffull[ChBF])[1], 10)
#     if pearson == 0.0:
#         plt.annotate(f"p < 0.0001 and r = {np.round(correlation, 5)}", xy=(250, 5),
#                      xycoords='axes points', fontsize=36, color='red')
#     else:
#         plt.annotate(f"p = {np.round(pearson, 5)} and r = {np.round(correlation, 5)}",
#                      xy=(250, 5), xycoords='axes points', fontsize=32, color='red')
#     #plt.title(f"PiB {study_cohort} {roi} and {cog} vol correlation", fontsize=20)
#     plt.savefig(f"{study_cohort}, Amyloid Association with {ChBF} Volume.png")
#     plt.clf()
#     plt.close(fig)
    
# for ChBF in ChBFlist:
#     fig, ax = plt.subplots(figsize=(14, 10), facecolor='white')
#     ax.scatter(x=cleandffull['age_at_visit'], y=cleandffull[ChBF], s=75)
#     ax.grid(axis="both")
#     ax.set_facecolor('white')
#     ax.spines['bottom'].set_color('k')
#     ax.spines['left'].set_color('k')
#     ax.spines['top'].set_color('k')
#     ax.spines['right'].set_color('k')
#     plt.xticks(fontsize=28)
#     plt.yticks(fontsize=28)
#     plt.xlabel("Global Amyloid (Centiloids)", fontsize=36)
#     plt.ylabel(f"{ChBF} Volume", fontsize=32)

#     # fit with polyfit
#     coef = np.polyfit(cleandffull['age_at_visit'], cleandffull[ChBF], 1)
#     poly1d_fn = np.poly1d(coef)
#     plt.plot(cleandffull['age_at_visit'], cleandffull[ChBF], 'bo', cleandffull['age_at_visit'], poly1d_fn(cleandffull['age_at_visit']), '-k')

#     correlation = np.corrcoef(cleandffull['age_at_visit'], cleandffull[ChBF])[0, 1]
#     pearson = np.round(stats.pearsonr(cleandffull['age_at_visit'], cleandffull[ChBF])[1], 10)
#     if pearson == 0.0:
#         plt.annotate(f"p < 0.0001 and r = {np.round(correlation, 5)}", xy=(250, 5),
#                      xycoords='axes points', fontsize=36, color='red')
#     else:
#         plt.annotate(f"p = {np.round(pearson, 5)} and r = {np.round(correlation, 5)}",
#                      xy=(250, 5), xycoords='axes points', fontsize=32, color='red')
#     #plt.title(f"PiB {study_cohort} {roi} and {cog} vol correlation", fontsize=20)
#     plt.savefig(f"{study_cohort}, Age Association with {ChBF} Volume.png")
#     plt.clf()
#     plt.close(fig)

cleandffull["de_gender"] = cleandffull["de_gender"].replace([1], 'Male')
cleandffull["de_gender"] = cleandffull["de_gender"].replace([2], 'Female')
cleandffull = cleandffull.rename({"de_gender":"Gender", "age_at_visit":"Age", "Left_ChBF": "Left ChBF Volume normalized to ICV", "Right_ChBF": "Right ChBF Volume normalized to ICV", "centiloid_value": "Amyloid (centiloids)", "tracer": "Tracer"}, axis=1)

ChBFlist = ('Left ChBF Volume normalized to ICV', 'Right ChBF Volume normalized to ICV')

for ChBF in ChBFlist:
    fig = sns.lmplot(data=cleandffull, x="Amyloid (centiloids)", y=ChBF)
    slope, intercept, r_value, p_value, std_err = stats.linregress(cleandffull['Amyloid (centiloids)'],cleandffull[ChBF])
    if p_value < 0.0001:
        plt.annotate(f"p < 0.0001 and r = {np.round(r_value, 5)}", xy=(5, 5),
                     xycoords='axes points', fontsize=8, color='red')
    else:
        plt.annotate(f"p = {np.round(p_value, 5)} and r = {np.round(r_value, 5)}",
                     xy=(5, 5), xycoords='axes points', fontsize=8, color='red')
    plt.savefig(f"{study_cohort}, Amyloid Association with {ChBF} Volume.png")
    plt.clf()
 
        
for ChBF in ChBFlist:
    fig = sns.lmplot(data=cleandffull, x="Age", y=ChBF)
    slope, intercept, r_value, p_value, std_err = stats.linregress(cleandffull['Age'],cleandffull[ChBF])
    if p_value < 0.0001:
        plt.annotate(f"p < 0.0001 and r = {np.round(r_value, 5)}", xy=(5, 5),
                     xycoords='axes points', fontsize=8, color='red')
    else:
        plt.annotate(f"p = {np.round(p_value, 5)} and r = {np.round(r_value, 5)}",
                     xy=(5, 5), xycoords='axes points', fontsize=8, color='red')
    plt.savefig(f"{study_cohort}, Age with {ChBF} Volume.png")
    plt.clf()

males_df = cleandffull.loc[cleandffull.Gender=="Male"]
females_df = cleandffull.loc[cleandffull.Gender=="Female"]

for ChBF in ChBFlist:
   fig = sns.lmplot(data=cleandffull, x="Amyloid (centiloids)", y=ChBF, hue="Gender")
   slope, intercept, r_value, p_value, std_err = stats.linregress(males_df['Amyloid (centiloids)'], males_df[ChBF])

   if p_value < 0.0001:
       plt.annotate(f"Males: p < 0.0001 and r = {np.round(r_value, 5)}", xy=(5, 15),
                 xycoords='axes points', fontsize=8, color='red')
   else:
       plt.annotate(f"Males: p = {np.round(p_value, 5)} and r = {np.round(r_value, 5)}",
                 xy=(5, 15), xycoords='axes points', fontsize=8, color='red')
   slope, intercept, r_value, p_value, std_err = stats.linregress(females_df['Amyloid (centiloids)'], females_df[ChBF])
  
   if p_value < 0.0001:
       plt.annotate(f"Females: p < 0.0001 and r = {np.round(r_value, 5)}", xy=(5, 5),
                     xycoords='axes points', fontsize=8, color='red')
   else:
       plt.annotate(f"Females: p = {np.round(p_value, 5)} and r = {np.round(r_value, 5)}",
                     xy=(5, 5), xycoords='axes points', fontsize=8, color='red')
   plt.savefig(f"{study_cohort}, Amyloid Association with {ChBF} Volume Separated by sex.png")
   plt.clf()
   
for ChBF in ChBFlist:

    fig = sns.relplot(data=cleandffull, x="Amyloid (centiloids)", y=ChBF, hue="Age", palette='Reds')
    sm = plt.cm.ScalarMappable(cmap="Reds")
    plt.savefig(f"{study_cohort}, Amyloid Association with {ChBF} Volume with age.png")
    plt.clf()
    
for ChBF in ChBFlist:
    fig = sns.lmplot(data=cleandffull, x="Age", y=ChBF, hue="Gender")
    slope, intercept, r_value, p_value, std_err = stats.linregress(males_df['Age'], males_df[ChBF])
    if p_value < 0.0001:
        plt.annotate(f"Males: p < 0.0001 and r = {np.round(r_value, 5)}", xy=(5, 15),
              xycoords='axes points', fontsize=8, color='red')
    else:
        plt.annotate(f"Males: p = {np.round(p_value, 5)} and r = {np.round(r_value, 5)}",
                 xy=(5, 15), xycoords='axes points', fontsize=8, color='red')
    slope, intercept, r_value, p_value, std_err = stats.linregress(females_df['Age'], females_df[ChBF])
  
    if p_value < 0.0001:
        plt.annotate(f"Females: p < 0.0001 and r = {np.round(r_value, 5)}", xy=(5, 5),
                     xycoords='axes points', fontsize=8, color='red')
    else:
        plt.annotate(f"Females: p = {np.round(p_value, 5)} and r = {np.round(r_value, 5)}",
                     xy=(5, 5), xycoords='axes points', fontsize=8, color='red')
    plt.savefig(f"{study_cohort}, Age association with {ChBF} Volume Separated by sex.png")
    plt.clf()
    
#    PiB_df = cleandffull.loc[cleandffull.Tracer=="PiB"]
#    FBP_df = cleandffull.loc[cleandffull.Tracer=="FBP"]
    
#for ChBF in ChBFlist:
#   fig = sns.lmplot(data=cleandffull, x="Amyloid (centiloids)", y=ChBF, hue="Tracer")
#   correlation = np.corrcoef(PiB_df['Amyloid (centiloids)'], PiB_df[ChBF])[0, 1]
#   pearson = np.round(stats.pearsonr(PiB_df['Amyloid (centiloids)'], PiB_df[ChBF])[1], 10)
#   if pearson == 0.0:
#       plt.annotate(f"PiB: p < 0.0001 and r = {np.round(correlation, 5)}", xy=(180, 5),
#                 xycoords='axes points', fontsize=8, color='red')
#   else:
#       plt.annotate(f"PiB: p = {np.round(pearson, 5)} and r = {np.round(correlation, 5)}",
#                 xy=(180, 15), xycoords='axes points', fontsize=8, color='red')
#   correlation = np.corrcoef(FBP_df['Amyloid (centiloids)'], FBP_df[ChBF])[0, 1]
#   pearson = np.round(stats.pearsonr(FBP_df['Amyloid (centiloids)'], FBP_df[ChBF])[1], 10)
#   if pearson == 0.0:
#       plt.annotate(f"FBP: p < 0.0001 and r = {np.round(correlation, 5)}", xy=(180, 5),
#                     xycoords='axes points', fontsize=8, color='red')
#   else:
#        plt.annotate(f"FBP: p = {np.round(pearson, 5)} and r = {np.round(correlation, 5)}",
#                     xy=(180, 5), xycoords='axes points', fontsize=8, color='red')
#   plt.savefig(f"{study_cohort}, Amyloid Association with {ChBF} Volume Separated by Tracer.png")
#   plt.clf()

x = datetime.datetime.now()
cleandffull.to_excel(f"ABC-DS ScLimbic DF for R {x.year}{x.month}{x.day}.xlsx")
