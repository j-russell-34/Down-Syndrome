# code to generate correlations between ROI specific FEOBV uptake, and cognitive performance

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
import scipy.stats as stats
import statsmodels.api as sm
import seaborn as sns
import functions






study_cohort = "DSCHOL"
#feobv_csv = "DSCHOL_FEOBV.csv"
cog_csv = "DSCHOL_Cognition.csv"


# import csvs
#feobv_df = pd.read_csv(feobv_csv)
cog_df = pd.read_csv(cog_csv)

# get subjects who have had FEOBV scans
#subject_list = functions.feobvsubs(feobv_df)

# Sort dataframes to same subject order

#feobv_df = feobv_df.sort_values("SUBJECT")

# drop unwanted columns
# cog_df = cog_df.drop(["record_id"], axis=1)

# Change NaNs to 0
# cog_df = cog_df.fillna(0)
# cog_df = cog_df.replace(to_replace=['<4'], value=[4])

# Remove spaces from cognition df and select subjects that underwent FEOBV:
cleancog_df = cog_df.groupby(["id"]).max().reset_index()

#feobvsub_cog_df = functions.select_feobv_sub(subject_list, df=cleancog_df, col_id="id")

#finalcog_df = feobvsub_cog_df.drop(["id"], axis=1)

# remove unwanted columns from FEOBV
#clean_feobv_df = feobv_df.drop(["SUBJECT", "ASSR", "PROJECT", "DATE", "SESSION", "SESSTYPE", "SITE", "PROCTYPE",
 #                               "cblmwm_suvr", "cortwm_eroded_suvr", "cortwm_suvr", "supravwm_eroded_suvr",
  #                              "supravwm_suvr"], axis=1)

#finalcog_df.index = np.arange(0, len(finalcog_df))


# Generate matrix for r and p-values for cog
#functions.genmatrix(clean_feobv_df, finalcog_df, study_cohort)

#functions.generate_scatter_p(clean_feobv_df, finalcog_df, study_cohort)




