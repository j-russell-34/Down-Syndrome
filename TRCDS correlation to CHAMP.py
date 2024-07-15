# code to correlate recent data to published ALbin 2018 FEOBV SUVRs in healthy volunteers

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy.stats as stats
import statsmodels.api as sm

compare = 'FEOBVQA_v2_TRCDS.csv'
studycohort = "Down Syndrome Individuals"
# for testing, Comment out and comment in line above when complete
# compare = "FEOBVQA_v2.csv"
# studycohort = "CHAMP"

# import csvs
newdatacsv = pd.read_csv(compare, encoding='utf-8')
champcsv = pd.read_csv('FEOBVQA_V2_CHAMP.csv', encoding='utf-8')

# remove any ROIs not in the Albin paper
champcsv = champcsv.dropna(axis=1,how='all')


comparisondf = pd.DataFrame([], columns =['ROI', studycohort, 'Normal Individuals'])

#drop unwanted columns
champcsv=champcsv.iloc[:,8:]
newdatacsv=newdatacsv.iloc[:,8:]


champcsv = champcsv.drop(['cblmwm_suvr','compositegm_suvr','cortwm_eroded_suvr', 'cortwm_suvr', 'supravwm_eroded_suvr', 'supravwm_suvr', 'ventraldc_suvr'], axis=1)
newdatacsv = newdatacsv.drop(['cblmwm_suvr','compositegm_suvr','cortwm_eroded_suvr', 'cortwm_suvr', 'supravwm_eroded_suvr', 'supravwm_suvr', 'ventraldc_suvr'], axis=1)


# iterate through
for roi in champcsv.columns:
    roimean = newdatacsv.loc[:, roi].mean()
    champ = champcsv.loc[:, roi].mean()
    add = pd.DataFrame([[roi, roimean, champ]], columns = ['ROI', studycohort, 'Normal Individuals'])
    #generate df with "ROI" "data set for comparison" "Albin dataset"
    comparisondf = pd.concat([comparisondf, add])



# linear regression
model = sm.OLS(comparisondf[studycohort], sm.add_constant(comparisondf['Normal Individuals'])).fit()
residuals = model.resid

# identify data points with large residuals
#outliers = np.where(np.abs(residuals) > 2*np.std(residuals))[0]
#outlier_labels = []
#for i in outliers:
    #outlier_to_app = comparisondf.iloc[i, comparisondf.columns.get_loc('ROI')]
    #outlier_labels.append(outlier_to_app)

# generate scatter plot
fig, ax = plt.subplots(figsize=(12, 12))
ax.scatter(x=comparisondf['Normal Individuals'], y=comparisondf[studycohort], s=150)
plt.xlabel("FEOBV SUVRs from Normal Individuals", fontsize= 20)
plt.ylabel(f"FEOBV SUVRs from {studycohort}", fontsize=20)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
ax.grid(axis="both")
ax.grid(False)
ax.set_facecolor('white')
ax.spines['bottom'].set_color('k')
ax.spines['left'].set_color('k')
ax.spines['top'].set_color('k')
ax.spines['right'].set_color('k')
#fit with polyfit
plt.plot(comparisondf['Normal Individuals'], model.predict(), color='k')


# highlight outliers in red and add labels
#for i, j in enumerate(outliers):
    #ax.scatter(comparisondf.iloc[j, comparisondf.columns.get_loc('albin2018')], comparisondf.iloc[j, comparisondf.columns.get_loc(studycohort)], color='r')
    #ax.annotate(comparisondf.iloc[j, comparisondf.columns.get_loc('ROI')], (comparisondf.iloc[j, comparisondf.columns.get_loc('albin2018')], comparisondf.iloc[j, comparisondf.columns.get_loc(studycohort)]))


# calculate correlations
correlation = np.corrcoef(comparisondf['Normal Individuals'], comparisondf[studycohort])[0,1]

pearson = np.round(stats.pearsonr(comparisondf['Normal Individuals'], comparisondf[studycohort])[1], 10)
if pearson == 0.0:
    plt.annotate(f"p < 0.0001 and r = {np.round(correlation, 5)}", fontsize=22, color='red', xy=(320, 20),
                 xycoords='axes points')
else:
    plt.annotate(f"p = {np.round(pearson,5)} and r = {np.round(correlation, 5)}", fontsize=22, color='red', xy=(320, 20),
                 xycoords='axes points')

#plt.title(f"{studycohort} comparison to Albin 2018", fontsize=20)
plt.savefig(f"{studycohort} comparison to Normal Individuals")

