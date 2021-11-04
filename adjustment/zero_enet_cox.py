from sksurv.linear_model import CoxnetSurvivalAnalysis
from sklearn.model_selection import KFold
from sksurv.metrics import concordance_index_censored
from scipy.optimize import minimize
from scipy.optimize import Bounds
import numpy as np
import pandas as pd
import pdb

y = pd.read_csv("y_train.csv")
x = pd.read_csv("train_df.gz")

x = x.to_numpy()
y = y.to_numpy()

#pdb.set_trace()

kf = KFold(n_splits=5,random_state=1,shuffle=True).split(y[:,1])

all_conc = []
all_coefs = []
all_alpha = []
n_fold = 0

for train_index, test_index in kf:
  print(n_fold)

  #pdb.set_trace()
  train_y = y[train_index,:]
  train_x = x[train_index,:]

  test_y = y[test_index,:]
  test_x = x[test_index,:]

  str_train_y = np.core.records.fromarrays(train_y[:,[1,0]].transpose(), names='time, status', formats = '?, i8')
  str_test_y = np.core.records.fromarrays(test_y[:,[1,0]].transpose(), names='time, status', formats = '?, i8')

  custom_alphas = np.logspace(1, 5, 50)/10000

  def custom_func(input_alpha):
    estimator = CoxnetSurvivalAnalysis(alphas = input_alpha, normalize=False, l1_ratio=0.5, verbose=True, fit_baseline_model=False, copy_X=False, max_iter = 1000, tol = 1e-04)

    estimator.fit(train_x, str_train_y)

    pred = estimator.predict(test_x, alpha = estimator.alphas_[0])

    curr_conc = concordance_index_censored(test_y[:,1].astype("bool"), test_y[:,0], pred)[0]

    return(curr_conc)


  start_tries = [0.0001, 0.001, 0.01, 0.1]
  def iterate_over(tries, last_conc, pass_through):
    print("start iterate over")
    #print(pass_through)
    back_conc = []
    for i in range(len(tries)):
      #print(tries[i])
      back_conc.append(custom_func(np.array([tries[i]])))
    #print("got all concs")
    back_conc = np.array(back_conc)
    max_ind = np.where(back_conc == max(back_conc))[0][0]
    #print("max_ind")
    #print(max_ind)
    if (abs(max(back_conc) - last_conc) < 0.001 and pass_through > 2) or pass_through > 5:
      return(tries[max_ind])
    else:
      if max_ind == (len(tries)-1):
        other_ind = max_ind-1
      elif max_ind == 0:
        other_ind = max_ind +1
      else:
        other_ind = np.array([max_ind-1,max_ind+1])[back_conc[[max_ind-1,max_ind+1]] == max(back_conc[[max_ind-1,max_ind+1]])][0]
      #print("other_ind")
      #print(other_ind)
      try_range = np.array(tries)[np.sort(np.array([max_ind, other_ind]))]
      new_tries = np.linspace(try_range[0], try_range[1], 5)
      #print(new_tries)
      return(iterate_over(tries, max(back_conc), pass_through+1))

  best_alpha = iterate_over(start_tries, 0, 1)

  #minimize(custom_func, np.array([0.01]), bounds = Bounds(0, 10), tol = 0.001)
  #array([0.00584834]

  estimator = CoxnetSurvivalAnalysis(alphas = [best_alpha], normalize=False, l1_ratio=0.5, verbose=True, fit_baseline_model=False, copy_X=False, max_iter = 1000, tol = 1e-04)
  estimator.fit(train_x, str_train_y)
  pred = estimator.predict(test_x, alpha = estimator.alphas_[0])
  curr_conc = concordance_index_censored(test_y[:,1].astype("bool"), test_y[:,0], pred)[0]

  all_conc.append(curr_conc)
  all_coefs.append(estimator.coef_)
  all_alpha.append(best_alpha)

  n_fold += 1


#pdb.set_trace()

#take weighted average
#all_conc = sum(all_conc)/len(all_conc)
#best_alpha = all_conc[all_conc[:,1] == max(all_conc[:,1]),0][0]
best_alpha = sum(np.array(all_conc)/sum(all_conc) * np.array(best_alpha))

#best_ind = np.where(all_conc[:,1] == max(all_conc[:,1]))[0][0]
#best_coefs = np.zeros((all_coefs[0].shape[0]))
#for i in range(len(all_coefs)):
#  best_coefs += all_coefs[i][:,best_ind]
#best_coefs = best_coefs/len(all_coefs)
#pd.DataFrame(best_coefs).to_csv("best_coef.csv")


estimator = CoxnetSurvivalAnalysis( normalize=False, l1_ratio=0.5, verbose=True, fit_baseline_model=True, copy_X=False, alphas = [best_alpha])
str_y = np.core.records.fromarrays(y[:,[1,0]].transpose(), names='time, status', formats = '?, i8')
estimator.fit(x, str_y)

pd.DataFrame(estimator.predict(x)).to_csv("train_pred.txt.gz", header = False)

del x
del y
y = pd.read_csv("y_test.csv")
x = pd.read_csv("test_df.gz")

static_pred = estimator.predict(x)

pd.DataFrame(estimator.coef_).to_csv("coef.txt.gz", header = False)
pd.DataFrame(static_pred).to_csv("static_pred.txt.gz", header = False)



#pred = estimator.predict_cumulative_hazard_function(x)
#del x

#input_time = np.arange(10, max(y["time"].to_numpy()), 10)
#cum_haz_pred = np.zeros((y.shape[0], len(input_time)))
#for i in range(cum_haz_pred.shape[0]):
#  cum_haz_pred[i,:] = pred[i](input_time)

#pd.DataFrame(cum_haz_pred).to_csv("cum_haz_pred.txt.gz", header = False)
#pd.DataFrame(input_time).to_csv("input_time.txt.gz", header = False)

print("done")
