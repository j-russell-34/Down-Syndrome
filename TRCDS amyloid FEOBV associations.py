import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
import scipy.stats as stats
import statsmodels.api as sm
import seaborn as sns
import functions

# import Amyvid and FEOBV uptakes

study_cohort = "DS_CHOL"
amyvid_csv = "DSCHOL_PiB.csv"
feobv_csv = "DSCHOL_FEOBV.csv"

radiotracer = "PiB"

amyvid_df = pd.read_csv(amyvid_csv)
feobv_df = pd.read_csv(feobv_csv)

# get subjects who had an amyvid scan
subject_list = functions.feobvsubs(amyvid_df, col_id="SUBJECT")

# get subjects who had an FEOBV scan
subject_list_feobv = functions.feobvsubs(feobv_df, col_id="SUBJECT")

# get just subjects with amyvid scans from feobv df, and deobv scans from amyvid df:

amyvidsub_feobv_df = functions.select_feobv_sub(subject_list=subject_list, df=feobv_df, col_id="SUBJECT")
feobvsub_amyvid_df = functions.select_feobv_sub(subject_list=subject_list_feobv, df=amyvid_df, col_id="SUBJECT")

# sort FEOBV and AMYVID dataframes to ensure order
amyvidsub_feobv_df = amyvidsub_feobv_df.sort_values("SUBJECT")
feobvsub_amyvid_df = feobvsub_amyvid_df.sort_values("SUBJECT")

# remove unwanted columns from FEOBV
clean_feobv_df = amyvidsub_feobv_df.drop(["SUBJECT", "ASSR", "PROJECT", "DATE", "SESSION", "SESSTYPE", "SITE", "PROCTYPE",
                                "cblmwm_suvr", "cortwm_eroded_suvr", "cortwm_suvr", "supravwm_eroded_suvr",
                                "supravwm_suvr"], axis=1)


# remove unwanted columns from AMYVID
clean_amyvid_df = feobvsub_amyvid_df.drop(["SUBJECT", "ASSR", "PROJECT", "DATE", "SESSION", "SESSTYPE", "SITE", "PROCTYPE",
                                "cblmwm_suvr", "cortwm_suvr"], axis=1)

clean_feobv_df.columns = clean_feobv_df.columns.str.strip("ROI_")
clean_feobv_df.columns = clean_feobv_df.columns.str.strip("_suvr")

clean_amyvid_df.columns = clean_amyvid_df.columns.str.strip("_suvr")

correlation_table=[]

clean_feobv_df = clean_feobv_df.rename({"hippocamp":"Hippocampus"}, axis=1)
clean_amyvid_df = clean_amyvid_df.rename({"hippocamp":"Hippocampus"}, axis=1)

for col in clean_feobv_df:
    correlation_table.append((col, "%.4f" % np.round((pearsonr(clean_amyvid_df[col], clean_feobv_df[col])[0]), 4),
                              "%.4f" % np.round((pearsonr(clean_amyvid_df[col], clean_feobv_df[col])[1]), 4)))

    if functions.pearsonr_pval(clean_amyvid_df[col], clean_feobv_df[col]) < 0.05:
        fig, ax = plt.subplots(figsize=(14, 10), facecolor='white')
        ax.scatter(x=clean_amyvid_df[col], y=clean_feobv_df[col], s=175, edgecolors=None, c="blue")
        ax.set_facecolor('white')
        ax.spines['bottom'].set_color('k')
        ax.spines['left'].set_color('k')
        ax.spines['top'].set_color('k')
        ax.spines['right'].set_color('k')
        plt.xticks(fontsize=28)
        plt.yticks(fontsize=28)
        plt.grid(False)
        plt.xlabel(f"{col} {radiotracer} SUVR", fontsize=36)
        plt.ylabel(f"{col} FEOBV SUVR", fontsize=32)

        # fit with polyfit
        coef = np.polyfit(clean_amyvid_df[col], clean_feobv_df[col], 1)
        poly1d_fn = np.poly1d(coef)
        plt.plot(clean_amyvid_df[col], clean_feobv_df[col], 'bo', clean_amyvid_df[col], poly1d_fn(clean_amyvid_df[col]),
                  '-k')

        correlation = np.corrcoef(clean_amyvid_df[col], clean_feobv_df[col])[0, 1]
        pearson = np.round(stats.pearsonr(clean_amyvid_df[col], clean_feobv_df[col])[1], 10)
        if pearson == 0.0:
            plt.annotate(f"p < 0.0001 and r = {np.round(correlation, 5)}", xy=(25, 25),
                          xycoords='axes points', fontsize=32, color='red')
        else:
            plt.annotate(f"p = {np.round(pearson, 5)} and r = {np.round(correlation, 5)}",
                          xy=(25, 25), xycoords='axes points', fontsize=32, color='red')
        #plt.title(f"{study_cohort} {col} FEOBV and {col} {radiotracer} SUVRs correlation", fontsize=18)
        plt.savefig(f"{study_cohort} {col} FEOBV and {col} {radiotracer} SUVRs correlation.png")
        plt.clf()
        plt.close(fig)

correlations = pd.DataFrame(correlation_table, columns=['ROI', 'Pearsons R', "p-value"])

functions.save_df_as_image_white(correlations, f"FEOBV and {radiotracer} SUVR relationships in {study_cohort}")

