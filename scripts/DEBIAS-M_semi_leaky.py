import sys
sys.path.append('/users/aflemis1/my_python_libs')

from debiasm import DebiasMClassifier
import debiasm
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

np.set_printoptions(threshold=sys.maxsize)

# parse folder argument
input = str(sys.argv[2])
output = str(sys.argv[3])

log_DEBIAS = pd.read_csv(f"./csv_files/{input}/lognorm_data.csv", dtype={'sample_name': str, 'Study_ID': str}) # read in data
data = log_DEBIAS.iloc[:,3:log_DEBIAS.shape[1]] # select the columns case, ID, and data columns (0 through amount of studies - 1)


## the read count matrix and "batches" a.k.a sample_IDs
X = data.iloc[:,2:data.shape[1]]
X = X.to_numpy()


## select the labels
y =  data["case"].to_numpy() # the 0s and 1s... make sure true (case/1) false(control/0)

# batches test
# batches = data.iloc[:,1:2]
# batches = batches.to_numpy()

# batches
batches = data.iloc[:, 1]
batches  = batches.to_numpy().astype(int)

## we assume the batches (Study_IDs) are numbered ints starting at '0',
## and they are in the first column of the input X matrices

# run this one without editing
# X_with_batch = np.hstack((batches[:, np.newaxis], X))
X_with_batch = np.hstack((batches[:, np.newaxis], X))
# X_with_batch[:10, :10]


## each time you're training without the study passed as a command line argument. If you do a LOO you'll have an amount of spreadsheets equal to the number of study_IDs * 2 
# (1 weight file and one "normalized" dataframe)

# this means the training will be done without the study you passed as a command line argument 
i = int(sys.argv[1]) 
i = i-1
val_inds = batches==i
# print(val_inds)
X_train, X_val = X_with_batch[~val_inds], X_with_batch[val_inds]
y_train, y_val = y[~val_inds], y[val_inds]


np.random.seed(123)
dmc = DebiasMClassifier(x_val=X_val) ## give it the held-out inputs to account for
                                    ## those domains shifts while training
np.random.seed(123)
dmc.fit(X_train, y_train)
print('finished training!')

print(roc_auc_score(y_val, dmc.predict_proba(X_val)[:, 1]))

X_debiassed = dmc.transform(X_with_batch)


#### write to CSV ####
id_dict = {}
Study_ID = log_DEBIAS["Study_ID"].unique()
id_array = log_DEBIAS["ID"].unique()
for id in range(0,len(Study_ID)):
    id_dict[id_array[id]] = Study_ID[id]

# get all the columns with only data (a.k.a bacteria)
cols = list(data.columns[2:data.shape[1]])


debiased_df = pd.DataFrame(data = X_debiassed,
                           columns = cols)


# merge the metadata columns with the data columns
first_cols = log_DEBIAS.iloc[:, 0:3]

frames = [first_cols, debiased_df]

merged = pd.concat(frames, axis=1)

merged.to_csv(f"./output/{output}/DEBIAS-M_runs/builtenv_debiased_lognorm_{id_dict[i]}.csv", index = False)


# get all the columns with only data (a.k.a bacteria)
weights = dmc.model.batch_weights
np_weights = weights.detach().numpy()
cols = list(data.columns[2:data.shape[1]])


# Study_ID = log_DEBIAS["Study_ID"].unique()

weights_df = pd.DataFrame(data = np_weights,
                          columns = cols,
                          index = Study_ID)

weights_df = weights_df.transpose()



weights_df.to_csv(f"./output/{output}/DEBIAS-M_runs/builtenv_debias_weights_{id_dict[i]}.csv")

print("done writing all files to csv")
