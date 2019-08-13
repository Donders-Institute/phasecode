import numpy as np
import scipy.io as sio
import scipy as scipy
import sklearn
from sklearn.linear_model import LogisticRegressionCV

project_dir = '/project/3011085.02/phasecode/';
results_dir = 'results/'
data = sio.loadmat(project_dir + results_dir + 'sub04_data_congruent.mat');

data1 = data['dat1'];
data2 = data['dat2'];
cov = data['cov']

ntrials = np.size(data1, 0)
design = np.concatenate((np.zeros((ntrials - 1, 1)), np.ones((ntrials - 1, 1))), axis=0);


tidx = 80;

prob = np.zeros((ntrials,2));
shortdata1 = data1[:,:,tidx]
#shortdata1 = np.reshape(shortdata1, (360, 41*300))
shortdata2 = data2[:,:,tidx]
#shortdata2 = np.reshape(shortdata2, (360, 41*300))

for trl in range(ntrials):
    print(trl)
    tmpdata1 = np.delete(shortdata1, trl, 0);
    tmpdata2 = np.delete(shortdata2, trl, 0);
    traindata = np.concatenate((tmpdata1, tmpdata2), axis=0);
    testdata = np.concatenate((np.reshape(shortdata1[trl, :], (1, -1)), np.reshape(shortdata2[trl, :], (1, -1))), axis=0)

    # prewhiten
    sigma = np.mean(np.delete(cov, trl, 0), 0);
    sigma_inv = scipy.linalg.fractional_matrix_power(sigma, -0.5)
    traindata = traindata @ sigma_inv
    testdata = testdata @ sigma_inv


    #scaler = sklearn.preprocessing.StandardScaler().fit(traindata)
    #traindata = scaler.transform(traindata)
    #test1 = scaler.transform(np.reshape(shortdata1[trl, :], (1, -1)))
    #test2 = scaler.transform(np.reshape(shortdata2[trl, :], (1, -1)))
    model = LogisticRegressionCV(penalty='elasticnet', max_iter=10000, solver='saga', l1_ratios=[0.1])
    model.fit(traindata, design.ravel())
    tmpprob = model.predict_proba(testdata)
    prob[trl-1,:]= np.append(tmpprob[0,0], tmpprob[1,1]);


print(prob1.mean())
print(prob2.mean())