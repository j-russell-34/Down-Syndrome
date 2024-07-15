import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
import scipy.stats as stats
import statsmodels.api as sm
import seaborn as sns





def feobvsubs(df, col_id):
    subject_list = []
    for index, row in df.iterrows():
        subject = [row[col_id]][0]
        subject_list.append(subject)
    return subject_list

def pearsonr_pval(x, y):
    return pearsonr(x, y)[1]

def save_df_as_image(df, path, title, mask=None):
    # Make plot
    sns.set(rc={'figure.figsize': (28, 22)})
    plot = sns.heatmap(df, mask=mask, cmap='coolwarm', annot=True, fmt='.4f')
    plt.title(title, fontsize=28)
    fig = plot.get_figure()
    fig.savefig(path)
    plt.clf()

def save_df_as_image_white(df, path):
    # Make plot
    fig, ax = plt.subplots(figsize=(10, 15))
    fig.patch.set_visible(False)
    ax.axis('off')
    ax.axis('tight')
    table = ax.table(cellText=df.values, colLabels=df.columns, loc='center')
    table.set_fontsize(14)
    table.scale(1, 1.25)
    fig.tight_layout()
    fig.savefig(path)
    plt.clf()

# def save_df_as_image_df2img(df, filename, title):
#     # Make plot
#     fig = df2img.plot_dataframe(
#         df,
#         title=dict(
#             font_color="black",
#             font_family="Ariel",
#             font_size=16,
#             text=title,
#         ),
#         tbl_header=dict(
#             align="center",
#             fill_color="white",
#             font_color="black",
#             font_size=10,
#             line_color="darkslategray",
#         ),
#         tbl_cells=dict(
#             align="right",
#             line_color="darkslategray",
#         ),
#         row_fill_color=("#ffffff", "#d7d8d6"),
#         fig_size=(300, 160),
#     )

#     df2img.save_dataframe(fig=fig, filename=filename)


def select_feobv_sub(subject_list, df, col_id):
    feobvsub_df = pd.DataFrame()
    for sub in subject_list:
        if type(df.loc[1, col_id]) == str:
            sub = str(sub)
        selected_sub = df.loc[(df[col_id] == sub)]
        feobvsub_df = pd.concat([feobvsub_df, selected_sub])
    feobvsub_df = feobvsub_df.sort_values(col_id)
    return feobvsub_df


def genmatrix(df1, df2, study_cohort):
    correlmatrix = pd.concat([df1.reset_index(drop=True), df2.reset_index(drop=True)], axis=1).reset_index()
    # Generate matrix for r and p-values for cog
    matrix = correlmatrix.corr().round(4).filter(df2.columns). \
        filter(df1.columns, axis=0)

    matrixps = correlmatrix.corr(method=pearsonr_pval).round(4). \
        filter(df2.columns).filter(df1.columns, axis=0)

    # create a mask based on the p-value matrix
    mask = matrixps >= 0.1

    # apply the mask to the original matrix

    save_df_as_image(matrix, f'{study_cohort} PiB SUVR and BCF vol R-values.png',
                     f'{study_cohort} PiB SUVR and BCF vol R-values')
    save_df_as_image(matrixps, f'{study_cohort} PiB SUVR and BCF vol p-values.png',
                     f'{study_cohort} PiB SUVR and BCF vol p-values')
    save_df_as_image(matrix, f'{study_cohort} PiB SUVR and BCF vol significant R-values.png',
                     f'{study_cohort} PiB SUVR and BCF vol significant R-values', mask=mask)

def generate_scatter_p(df1,df2, study_cohort):
    for roi in df1:
        roi_namer = roi[:-5]
        if roi_namer == "postcing":
            roi_name = "Posterior Cingulate Cortex"
        else:
            roi_name = roi_namer[4:]
        for cog in df2:
            if cog == "dems_age":
                cog_name = "Age"
            else:
                chnuc = cog[2:-6]
                if cog[-5] == "L":
                    side = "Left"
                else:
                    side = "Right"
                cog_name = f"Basal Cholinergic Forebrain {side} Ch{chnuc} Volume"
            new_df2 = df2
            new_df1 = df1
            if df2[cog].isna().any():
                nanindex = list(zip(*np.where(np.asanyarray(np.isnan(df2[cog])))))
                for n in nanindex:
                    new_df2 = new_df2.drop([n[0]], axis=0)
                    new_df1 = new_df1.drop([n[0]], axis=0)
            if pearsonr_pval(new_df2[cog], new_df1[roi]) < 0.1:
                fig, ax = plt.subplots(figsize=(14, 10), facecolor='white')
                ax.scatter(x=new_df2[cog], y=new_df1[roi], s=175)
                #ax.grid(axis="false")
                ax.set_facecolor('white')
                ax.spines['bottom'].set_color('k')
                ax.spines['left'].set_color('k')
                ax.spines['top'].set_color('k')
                ax.spines['right'].set_color('k')
                #avoid overlapping - comment out if not needed
                #plt.setp(ax.get_xticklabels(), rotation=30, ha='right')
              
                plt.xticks(fontsize=20)
                plt.yticks(fontsize=28)
                plt.xlabel(f"{cog_name}", fontsize=32)
                plt.ylabel(f"{roi_name} FEOBV SUVR", fontsize=32)
                plt.tight_layout()
                # fit with polyfit
                coef = np.polyfit(new_df2[cog], new_df1[roi], 1)
                poly1d_fn = np.poly1d(coef)
                plt.plot(new_df2[cog], new_df1[roi], 'bo', new_df2[cog], poly1d_fn(new_df2[cog]), '-k')

                correlation = np.corrcoef(new_df2[cog], new_df1[roi])[0, 1]
                pearson = np.round(stats.pearsonr(new_df2[cog], new_df1[roi])[1], 10)
                if pearson == 0.0:
                    plt.annotate(f"p < 0.0001 and r = {np.round(correlation, 5)}", xy=(25, 25),
                                 xycoords='axes points', fontsize=36, color='red')
                else:
                    plt.annotate(f"p = {np.round(pearson, 5)} and r = {np.round(correlation, 5)}",
                                 xy=(25, 25), xycoords='axes points', fontsize=32, color='red')
                #plt.title(f"PiB {study_cohort} {roi} and {cog} vol correlation", fontsize=20)
                plt.savefig(f"FEOBV {study_cohort} {roi_name} SUVRs and {cog_name} correlation.png")
                plt.clf()
                plt.close(fig)


def correlpet(csvx, radiotracerx, csvy, radiotracery, study_cohort):
    # import Amyvid and FEOBV uptakes

    amyvid_csv = csvy
    feobv_csv = csvx

    radiotracer = radiotracery

    amyvid_df = pd.read_csv(amyvid_csv)
    feobv_df = pd.read_csv(feobv_csv)

    # get subjects who had an amyvid scan
    subject_list = feobvsubs(amyvid_df)

    # get subjects who had an FEOBV scan
    subject_list_feobv = feobvsubs(feobv_df)

    # get just subjects with amyvid scans from feobv df, and deobv scans from amyvid df:

    amyvidsub_feobv_df = select_feobv_sub(subject_list=subject_list, df=feobv_df, col_id="SUBJECT")
    feobvsub_amyvid_df = select_feobv_sub(subject_list=subject_list_feobv, df=amyvid_df, col_id="SUBJECT")

    # sort FEOBV and AMYVID dataframes to ensure order
    amyvidsub_feobv_df = amyvidsub_feobv_df.sort_values("SUBJECT")
    feobvsub_amyvid_df = feobvsub_amyvid_df.sort_values("SUBJECT")

    # remove unwanted columns from FEOBV
    clean_feobv_df = amyvidsub_feobv_df.drop(
        ["SUBJECT", "ASSR", "PROJECT", "DATE", "SESSION", "SESSTYPE", "SITE", "PROCTYPE",
         ], axis=1)

    for col in clean_feobv_df:
        if "wm_suvr" in col:
            clean_feobv_df = clean_feobv_df.drop([col], axis=1)
        if "eroded_suvr" in col:
            clean_feobv_df = clean_feobv_df.drop([col], axis=1)


    # remove unwanted columns from AMYVID
    clean_amyvid_df = feobvsub_amyvid_df.drop(
        ["SUBJECT", "ASSR", "PROJECT", "DATE", "SESSION", "SESSTYPE", "SITE", "PROCTYPE",
         ], axis=1)

    for hed in clean_amyvid_df:
        if "wm_suvr" in hed:
            clean_amyvid_df = clean_amyvid_df.drop([hed], axis=1)
        if "eroded_suvr" in col:
            clean_amyvid_df = clean_amyvid_df.drop([hed], axis=1)

    clean_feobv_df.columns = clean_feobv_df.columns.str.strip("ROI_")
    clean_feobv_df.columns = clean_feobv_df.columns.str.strip("_suvr")
    clean_amyvid_df.columns = clean_amyvid_df.columns.str.strip("_suvr")

    correlation_table = []

    for col in clean_feobv_df:
        correlation_table.append(
            (col, "%.4f" % np.round((pearsonr(clean_amyvid_df[col], clean_feobv_df[col])[0]), 4),
             "%.4f" % np.round((pearsonr(clean_amyvid_df[col], clean_feobv_df[col])[1]), 4)))

        if pearsonr_pval(clean_amyvid_df[col], clean_feobv_df[col]) < 0.05:
            fig, ax = plt.subplots(figsize=(10, 10))
            ax.scatter(x=clean_amyvid_df[col], y=clean_feobv_df[col])
            plt.xlabel(f"{col} {radiotracer} SUVR", fontsize=16)
            plt.ylabel(f"{col} {radiotracerx} SUVR", fontsize=16)

            # fit with polyfit
            coef = np.polyfit(clean_amyvid_df[col], clean_feobv_df[col], 1)
            poly1d_fn = np.poly1d(coef)
            plt.plot(clean_amyvid_df[col], clean_feobv_df[col], 'bo', clean_amyvid_df[col],
                     poly1d_fn(clean_amyvid_df[col]),
                     '-k')

            correlation = np.corrcoef(clean_amyvid_df[col], clean_feobv_df[col])[0, 1]
            pearson = np.round(stats.pearsonr(clean_amyvid_df[col], clean_feobv_df[col])[1], 10)
            if pearson == 0.0:
                plt.annotate(f"p < 0.0001 and r = {np.round(correlation, 5)}", xy=(400, 5),
                             xycoords='axes points', fontsize=12, color='red')
            else:
                plt.annotate(f"p = {np.round(pearson, 5)} and r = {np.round(correlation, 5)}",
                             xy=(400, 5), xycoords='axes points', fontsize=12, color='red')
            plt.title(f"{study_cohort} {col} {radiotracerx} and {col} {radiotracer} SUVRs correlation", fontsize=18)
            plt.savefig(f"{study_cohort} {col} {radiotracerx} and {col} {radiotracer} SUVRs correlation.png")
            plt.clf()
            plt.close(fig)

    correlations = pd.DataFrame(correlation_table, columns=['ROI', 'Pearsons R', "p-value"])

    save_df_as_image_white(correlations, f"{radiotracerx} and {radiotracer} SUVR relationships in {study_cohort}")




