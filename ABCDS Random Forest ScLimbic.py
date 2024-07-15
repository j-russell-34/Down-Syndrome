# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 13:25:25 2024

@author: russj13
"""

import pandas as pd
import numpy as np
import random
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score, confusion_matrix, precision_score, recall_score, ConfusionMatrixDisplay
from sklearn.model_selection import RandomizedSearchCV, train_test_split

import matplotlib.pyplot as plt
from sklearn.preprocessing import OneHotEncoder


sclimbic_df = pd.read_excel("sclimbic_abcds.xlsx")

#dnseg_df["Sex"] = dnseg_df["Sex"].map({'Female':1, 'Male':2})

X = sclimbic_df.drop(['SUBJECT','Gender','Tracer', 'ASSR.x', 'SESSION',
                  'PROJECT.x','PROCTYPE.x','totalBFavg', 'samseg_lesions'], axis=1)

X.rename(columns={'Amyloid..centiloids.': 'Amyloid Centiloids', 'samseg_sbtiv': 'ICV'}, inplace=True)

y = sclimbic_df["totalBFavg"]

onehotencoder = OneHotEncoder()

X1 = onehotencoder.fit_transform(X.Sex.values.reshape(-1,1)).toarray()

df = pd.DataFrame(X1, columns =['Female', 'Male']) 

X = pd.concat([X, df], axis=1)

X = X.drop(['Sex'], axis=1) 
 
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, 
                                                    random_state=3)

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
ax.set_title("Feature importances on BF Volume")
ax.set_ylabel("Importance score")
fig.tight_layout()
#plt.annotate("Variance Explained by Model: 13.9%",
 #             xy=(170, 200), xycoords='axes points', fontsize=8, color='red')
plt.savefig("RandomForest Total BF VarImp.png")
plt.clf()
plt.close(fig)
