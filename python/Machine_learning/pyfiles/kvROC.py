#!/usr/bin/env python
print(__doc__)

import pandas as pd
import numpy as np
from scipy import interp
import matplotlib.pyplot as plt
from sklearn import svm, datasets
from sklearn.metrics import auc
from sklearn.metrics import plot_roc_curve
from sklearn.model_selection import StratifiedKFold

iris = datasets.load_iris()


def kvroc(df, vr, ft, k):
    n_samples, n_features = df.shape
    random_state = np.random.RandomState(0)
    vr = np.c_[vr, random_state.randn(n_samples, 200 * n_features)]
    cv = StratifiedKFold(n_splits=k)
    classifier = svm.SVC(kernel='linear', probability=True,
                         random_state=random_state)
    tprs = []
    aucs = []
    mean_fpr = np.linspace(0, 1, 100)

    fig, ax = plt.subplots()
    for i, (train, test) in enumerate(cv.split(vr, ft)):
        classifier.fit(vr[train], ft[train])
        viz = plot_roc_curve(classifier, vr[test], ft[test],
                             name='ROC fold {}'.format(i),
                             alpha=0.3, lw=1, ax=ax)
        interp_tpr = interp(mean_fpr, viz.fpr, viz.tpr)
        interp_tpr[0] = 0.0
        tprs.append(interp_tpr)
        aucs.append(viz.roc_auc)

    ax.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r',
            label='Chance', alpha=.8)

    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(aucs)
    ax.plot(mean_fpr, mean_tpr, color='b',
            label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc),
            lw=2, alpha=.8)

    std_tpr = np.std(tprs, axis=0)
    tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
    tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
    ax.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
                    label=r'$\pm$ 1 std. dev.')

    ax.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05],
           title="Receiver operating characteristic example")
    ax.legend(loc="lower right")
    plt.show()
    fig.savefig("./results/kfoldcv.pdf")
    return fig


import os
os.listdir("./pythondata/")

data = pd.read_csv("./pythondata/mutation_model-multiCOX_risk_KIRC.csv")
x = data["status"]
y = data["risk"]
fig1 = kvroc(data, x, y, 10)

data = pd.read_csv("./pythondata/mutation_model-multiCOX_risk_GBM.csv")
x = data["status"]
y = data["risk"]
fig2 = kvroc(data, x, y, 10)

data = pd.read_csv("./pythondata/mutation_model-multiCOX_risk_DLBC.csv")
x = data["status"]
y = data["risk"]
fig3 = kvroc(data, x, y, 10)

data = pd.read_csv("./pythondata/mutation_model-multiCOX_risk_LUSC.csv")
x = data["status"]
y = data["risk"]
fig4 = kvroc(data, x, y, 10)

data = pd.read_csv("./pythondata/mutation_model-multiCOX_risk_ACC.csv")
x = data["status"]
y = data["risk"]
fig5 = kvroc(data, x, y, 10)

data = pd.read_csv("./pythondata/mutation_model-multiCOX_risk_BLCA.csv")
x = data["status"]
y = data["risk"]
fig6 = kvroc(data, x, y, 10)


data = pd.read_csv("./pythondata/mutation_model-multiCOX_risk_BRCA.csv")
x = data["status"]
y = data["risk"]
fig7 = kvroc(data, x, y, 10)


data = pd.read_csv("./pythondata/mutation_model-multiCOX_risk_CHOL.csv")
x = data["status"]
y = data["risk"]
fig8 = kvroc(data, x, y, 10)


data = pd.read_csv("./pythondata/mutation_model-multiCOX_risk_COAD.csv")
x = data["status"]
y = data["risk"]
fig9 = kvroc(data, x, y, 10)


data = pd.read_csv("./pythondata/mutation_model-multiCOX_risk_LAML.csv")
x = data["status"]
y = data["risk"]
fig10 = kvroc(data, x, y, 10)


data = pd.read_csv("./pythondata/mutation_model-multiCOX_risk_ESCA.csv")
x = data["status"]
y = data["risk"]
fig11 = kvroc(data, x, y, 10)


data = pd.read_csv("./pythondata/mutation_model-multiCOX_risk_KICH.csv")
x = data["status"]
y = data["risk"]
fig12 = kvroc(data, x, y, 10)


data = pd.read_csv("./pythondata/mutation_model-multiCOX_risk_KIRP.csv")
x = data["status"]
y = data["risk"]
fig13 = kvroc(data, x, y, 10)


data = pd.read_csv("./pythondata/mutation_model-multiCOX_risk_LUAD.csv")
x = data["status"]
y = data["risk"]
fig14 = kvroc(data, x, y, 10)


data = pd.read_csv("./pythondata/mutation_model-multiCOX_risk_MESO.csv")
x = data["status"]
y = data["risk"]
fig15 = kvroc(data, x, y, 10)


data = pd.read_csv("./pythondata/mutation_model-multiCOX_risk_OV.csv")
x = data["status"]
y = data["risk"]
fig16 = kvroc(data, x, y, 10)

data = pd.read_csv("./pythondata/mutation_model-multiCOX_risk_age_stage_MRscore.csv")
x = data["OS"]
y = data["risk"]
fig17 = kvroc(data, x, y, 10)

data = pd.read_csv("./pythondata/mutation_model-multiCOX_risk_filter.csv")
x = data["OS"]
y = data["risk"]
fig18 = kvroc(data, x, y, 10)