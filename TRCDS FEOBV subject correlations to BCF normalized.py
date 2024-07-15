# code to generate correlations between ROI specific FEOBV uptake, age and CBF vols

import pandas as pd
import functions



study_cohort = "Down Syndrome DSCHOL"
feobv_csv = "DSCHOL_FEOBV.csv"
demographics_csv = "DSChol_dems.csv"
cbf_csv = "DSCHOL_CBF.csv"
ticv_csv = "DSCHOL_samseg.csv"



# import csvs
feobv_df = pd.read_csv(feobv_csv)
demo_df = pd.read_csv(demographics_csv)
cbf_df = pd.read_csv(cbf_csv)
ticv_df = pd.read_csv(ticv_csv)

# check through cholinergic basal forebrain df and demographics df to obtain headers for correlation matrix
cbf_regions = []
for roi in cbf_df:
    if "CH" in roi and "VOL" in roi:
        cbf_regions.append(roi)


# # Generate df with column headings for correlation matrix
# correl_df = pd.DataFrame([], columns = ["ROI", "Age r-value", "Age-p-value"])
#
# for roi in cbf_regions:
#     correl_df[f"{roi} r-value"] = []
#     correl_df[f"{roi} p-value"] = []


# get subjects who had an FEOBV scan
subject_list = functions.feobvsubs(df=feobv_df, col_id="SUBJECT")

# get just subjects with FEOBV scans from BFC df and ticv df:

feobvsub_cbf_df = functions.select_feobv_sub(subject_list, df=cbf_df, col_id="SUBJECT")
feobvsub_ticv_df = functions.select_feobv_sub(subject_list, df=ticv_df, col_id="SUBJECT")

# sort FEOBV CHBF and TICV dataframes to ensure order
feobv_df = feobv_df.sort_values("SUBJECT")
feobvsub_cbf_df = feobvsub_cbf_df.sort_values("SUBJECT")
feobvsub_ticv_df = feobvsub_ticv_df.sort_values("SUBJECT")

# average CH BF vols where multiple scans:
feobvsub_cbf_df = feobvsub_cbf_df.drop(["ASSR", "PROJECT", "DATE", "SESSION", "SESSTYPE", "SITE", "PROCTYPE"], axis=1)
feobvsub_cbf_df = feobvsub_cbf_df[feobvsub_cbf_df.columns.drop(list(feobvsub_cbf_df.filter(regex='MEAN')))]
feobvsub_cbf_df = feobvsub_cbf_df[feobvsub_cbf_df.columns.drop(list(feobvsub_cbf_df.filter(regex='TOT')))]
averagedcbf_df = feobvsub_cbf_df.groupby(["SUBJECT"]).mean().reset_index()
averagedcbf_df = averagedcbf_df.drop(["SUBJECT"], axis=1)

# normalize CH BF vols to total ICV
normalized_cbf = averagedcbf_df[["CH123_L_VOL", "CH123_R_VOL", "CH4_L_VOL","CH4_R_VOL"]].div(feobvsub_ticv_df["samseg_sbtiv"],axis=0)
        


# remove unwanted columns from FEOBV
clean_feobv_df = feobv_df.drop(["SUBJECT", "ASSR", "PROJECT", "DATE", "SESSION", "SESSTYPE", "SITE", "PROCTYPE",
                                "cblmwm_suvr", "cortwm_eroded_suvr", "cortwm_suvr", "supravwm_eroded_suvr",
                                "supravwm_suvr"], axis=1)

# select only subjects that have had FEOBV scans (in subjects list from demographics CSV
feobvsub_demo_df = functions.select_feobv_sub(subject_list, df=demo_df, col_id="id")

#sort demographics to same order as FEOBV and BCF
feobvsub_demo_df = feobvsub_demo_df.sort_values("id")

# select ages and add to dataframe for correlation matrix
agecolumn = feobvsub_demo_df.filter(regex='age')
df_full = pd.concat([normalized_cbf.reset_index(drop=True), agecolumn.reset_index(drop=True)], axis=1)

functions.genmatrix(clean_feobv_df, df_full, study_cohort)


# correlate FEOBV uptake with CH BF Vols and plot if p<0.05

functions.generate_scatter_p(clean_feobv_df, df_full, study_cohort)





