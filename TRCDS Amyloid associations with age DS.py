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
amyvid_csv = "AMYVIDQA_v2_Trc_ds.csv"
demo_csv = "DSChol_dems.csv"

radiotracer = "PiB"

amyvid_df = pd.read_csv(amyvid_csv)
demo_df = pd.read_csv(demo_csv)

# get subjects who had an amyvid scan
subject_list = functions.feobvsubs(amyvid_df)

# get just subjects with amyvid scans from demo df, and deobv scans from amyvid df:

amyvidsub_demo_df = functions.select_feobv_sub(subject_list=subject_list, df=demo_df, col_id="id")

# sort FEOBV and demo dataframes to ensure order
amyvidsub_demo_df = amyvidsub_demo_df.sort_values("id")
amyvid_df = amyvid_df.sort_values("SUBJECT")

# remove unwanted columns from FEOBV
clean_demo_df = amyvidsub_demo_df.drop(["id"], axis=1)


# remove unwanted columns from AMYVID
clean_amyvid_df = amyvid_df.drop(["SUBJECT", "ASSR", "PROJECT", "DATE", "SESSION", "SESSTYPE", "SITE", "PROCTYPE",
                                "cblmwm_suvr", "cortwm_suvr"], axis=1)


correlation_table=[]

for col in clean_amyvid_df:
    correlation_table.append((col, "%.4f" % np.round((pearsonr(clean_demo_df["dems_age"], clean_amyvid_df[col])[0]), 4),
                             "%.4f" % np.round((pearsonr(clean_demo_df["dems_age"], clean_amyvid_df[col])[1]), 4)))

    if functions.pearsonr_pval(clean_demo_df["dems_age"], clean_amyvid_df[col]) < 0.05:
        fig, ax = plt.subplots(figsize=(10, 10))
        ax.scatter(x=clean_demo_df["dems_age"], y=clean_amyvid_df[col])
        plt.xlabel(f"Participant Age", fontsize=16)
        plt.ylabel(f"{col} {radiotracer} SUVR", fontsize=16)

        # fit with polyfit
        coef = np.polyfit(clean_demo_df["dems_age"], clean_amyvid_df[col], 1)
        poly1d_fn = np.poly1d(coef)
        plt.plot(clean_demo_df["dems_age"], clean_amyvid_df[col], 'bo', clean_demo_df["dems_age"],
                 poly1d_fn(clean_demo_df["dems_age"]),'-k')

        correlation = np.corrcoef(clean_demo_df["dems_age"], clean_amyvid_df[col])[0, 1]
        pearson = np.round(stats.pearsonr(clean_demo_df["dems_age"], clean_amyvid_df[col])[1], 10)
        if pearson == 0.0:
            plt.annotate(f"p < 0.0001 and r = {np.round(correlation, 5)}", xy=(400, 5),
                         xycoords='axes points', fontsize=12, color='red')
        else:
            plt.annotate(f"p = {np.round(pearson, 5)} and r = {np.round(correlation, 5)}",
                         xy=(400, 5), xycoords='axes points', fontsize=12, color='red')
        plt.title(f"{study_cohort} {col} {radiotracer} SUVRs and age correlation", fontsize=18)
        plt.savefig(f"{study_cohort} {col} {radiotracer} SUVRs and age correlation.png")
        plt.clf()
        plt.close(fig)

correlations = pd.DataFrame(correlation_table, columns=['ROI', 'Pearsons R', "p-value"])

functions.save_df_as_image_white(correlations, f"{radiotracer} SUVR and age relationships in {study_cohort}")

