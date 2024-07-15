# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import pandas as pd
import numpy as np
import random
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score, confusion_matrix, precision_score, recall_score, ConfusionMatrixDisplay
from sklearn.model_selection import RandomizedSearchCV, train_test_split
from scipy.stats import randint
from sklearn.tree import export_graphviz
from IPython.display import Image
import graphviz
import matplotlib.pyplot as plt
from sklearn.preprocessing import OneHotEncoder



dnseg_df = pd.read_excel("dnseg_abcds.xlsx")

#dnseg_df["Sex"] = dnseg_df["Sex"].map({'Female':1, 'Male':2})

X = dnseg_df.drop(['SUBJECT','Gender','Tracer','Left_DnSeg', 'Right_DnSeg', 'ASSR',
                  'PROJECT','PROCTYPE','totalch4avg', 'avgch4divtiv', 'samseg_lesions'], axis=1)

X.rename(columns={'Amyloid..centiloids.': 'Amyloid Centiloids', 'samseg_sbtiv': 'ICV'}, inplace=True)

y = dnseg_df["totalch4avg"]

onehotencoder = OneHotEncoder()

X1 = onehotencoder.fit_transform(X.Sex.values.reshape(-1,1)).toarray()

df = pd.DataFrame(X1, columns =['Female', 'Male']) 

X = pd.concat([X, df], axis=1)

X = X.drop(['Sex'], axis=1) 
 
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=3)

rf = RandomForestRegressor(n_estimators=5000, criterion="squared_error", 
                           max_depth=7, min_samples_split=5, min_samples_leaf=2,
                           random_state=3)
rf.fit(X_train, y_train)

y_pred = rf.predict(X_test)

feature_names = ['Age', 'Amyloid Centiloids', 'ICV', 'Female', 'Male']

accuracy = r2_score(y_test, y_pred)
print("R2 of model:", accuracy)
importances = rf.feature_importances_
std = np.std([rf.feature_importances_ for tree in rf.estimators_], axis=0)
forest_importances = pd.Series(importances, index=feature_names)

fig, ax = plt.subplots()
forest_importances.plot.bar(yerr=std, ax=ax)
ax.set_title("Feature importances on Ch4 Volume")
ax.set_ylabel("Importance score")
fig.tight_layout()
#plt.annotate("Variance Explained by Model: 13.9%",
             # xy=(170, 200), xycoords='axes points', fontsize=8, color='red')
plt.savefig("RandomForest DnSeg VarImp.png")
plt.clf()
plt.close(fig)


