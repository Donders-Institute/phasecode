import numpy as np
import scipy.io as sio
import sklearn
from sklearn.linear_model import LogisticRegressionCV

project_dir = '/project/3011085.02/phasecode/';
results_dir = 'results/'
data = sio.loadmat(project_dir + results_dir + 'test4.mat');

data1 = data['t1'];
data2 = data['t2'];

design = np.concatenate((np.zeros((np.size(data1, 0) - 1, 1)), np.ones((np.size(data1, 0) - 1, 1))), axis=0);
ntrials = np.size(data1, 0)

tidx = 80;

prob1 = [];
prob2 = [];
shortdata1 = data1[:,:,tidx]
#shortdata1 = np.reshape(shortdata1, (360, 41*300))
shortdata2 = data2[:,:,tidx]
#shortdata2 = np.reshape(shortdata2, (360, 41*300))

for trl in range(ntrials):
    print(trl)
    tmpdata1 = np.delete(shortdata1, trl, 0);
    tmpdata2 = np.delete(shortdata2, trl, 0);
    traindata = np.concatenate((tmpdata1, tmpdata2), axis=0);
    scaler = sklearn.preprocessing.StandardScaler().fit(traindata)
    traindata = scaler.transform(traindata)
    test1 = scaler.transform(np.reshape(shortdata1[trl, :], (1, -1)))
    test2 = scaler.transform(np.reshape(shortdata2[trl, :], (1, -1)))
    model = LogisticRegressionCV(penalty='elasticnet', cv=5, max_iter=1000, solver='saga', l1_ratios=[0.1])
    model.fit(traindata, design)
    tmpprob1 = model.predict_proba(test1)
    tmpprob1 = tmpprob1[0];
    prob1 = np.append(prob1, tmpprob1[0])
    tmpprob2 = model.predict_proba(test2);
    tmpprob2 = tmpprob2[0];
    prob2 = np.append(prob2, tmpprob2[1])

print(prob1.mean())
print(prob2.mean())