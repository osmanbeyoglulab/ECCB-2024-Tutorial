import numpy as np
import pandas as pd
import scanpy as sc
import random, time
from scipy.stats import spearmanr, pearsonr


def assign_folds(adata, n_folds=5, train_percent=None, random_seed=0):
    random.seed(random_seed)
    if train_percent is not None:
        adata.var["fold"]=np.random.uniform(0,1, adata.n_vars)<train_percent
    else:
        adata.var["fold"]=np.random.randint(0,n_folds, adata.n_vars)
        

class ModelBase:
    def __init__(self, adata, layer="dca", gene_tf=None, intercept=False):
        self.adata = adata
        self.intercept = intercept
        if gene_tf is None:
            gene_tf = adata.varm['gene_tf']

        self.D = gene_tf
        self.Y = adata.to_df(layer).T
        self.genes = np.intersect1d(self.Y.index, self.D.index)
        self.D = self.D.loc[self.genes]
        self.Y = self.Y.loc[self.genes]

        self.n_folds = np.max(adata.var['fold']) + 1
        self.n_genes = len(self.genes)
        self.n_spots = adata.n_obs
        self.n_tfs = self.D.shape[1]

        self.training_genes = [adata.var.query("fold != @i").index for i in range(self.n_folds)]
        self.testing_genes = [adata.var.query("fold == @i").index for i in range(self.n_folds)]
        self.D_train = [self.D.loc[self.training_genes[i]].to_numpy() for i in range(self.n_folds)]
        self.D_test = [self.D.loc[self.testing_genes[i]].to_numpy() for i in range(self.n_folds)]
        self.Y_train = [self.Y.loc[self.training_genes[i]].to_numpy() for i in range(self.n_folds)]
        self.Y_test = [self.Y.loc[self.testing_genes[i]].to_numpy() for i in range(self.n_folds)]
        self.D = self.D.to_numpy()

        self.K = adata.obsm.get("kernel")
        self.svdD = []
        for i in range(self.n_folds):
            self.svdD.append(np.linalg.svd(self.D_train[i], full_matrices=False))
        self.svdK = None

        self.W = [None] * self.n_folds
        self.y_pred = [None] * self.n_folds


    def evaluate(self, fold=0, gene_set="testing"):
        if fold == -1:
            y_pred = np.asarray(self.y_pred_concat)
            Y = np.asarray(self.Y)
        else:
            W = self.W[fold]
            if gene_set == "testing":
                Y = self.Y_test[fold]
                D = self.D_test[fold]
            elif gene_set == "training":
                Y = self.Y_train[fold]
                D = self.D_train[fold]
            y_pred = D.dot(W)

        cor = [pearsonr(y_pred[:, spot], Y[:, spot])[0] for spot in range(self.n_spots)]
        gene_cor = [pearsonr(y_pred[i, :].T, Y[i, :].T)[0] for i in range(Y.shape[0])]
        return cor, gene_cor


    def fit(self, grid_search_params=None, fixed_params=dict(), verbose=False, stages=5, n_steps=5, axis=0):
        t1 = time.time()
        if grid_search_params is None:
            self.params = fixed_params
        else:
            self.params = self.grid_search(n_steps, grid_search_params, fixed_params=fixed_params, verbose=verbose, stages=stages, axis=axis)
        
        params = self.params
        for fold in range(self.n_folds):
            self.train(fold=fold, **params)
            self.y_pred[fold] = self.D_test[fold].dot(self.W[fold])

        self.W_concat = np.mean(self.W, axis=0)
        self.y_pred_concat = pd.DataFrame(data=None, index=self.Y.index, columns=self.Y.columns)

        for fold in range(self.n_folds):
            self.y_pred_concat.loc[self.testing_genes[fold]] = self.y_pred[fold]
        t2 = time.time()
        print('Time elapsed: %.2f seconds'%(t2-t1))


    def grid_search(self, n_steps, params, fixed_params=dict(), verbose=False, stages=2, axis=0):
        if verbose:
            print("stages remaining: " + str(stages))
        param_names = list(params.keys())
        n_params = len(param_names)
        param_steps = [np.log10(params[param_names[i]][1] / params[param_names[i]][0]) / (n_steps - 1) for i in
                       range(n_params)]
        param_mins = [np.log10(params[param_names[i]][0]) for i in range(n_params)]
        perf = list()
        best_perf = -1
        best_params= None
        if verbose:
            pbar = tqdm(range(n_steps ** n_params))
        else:
            pbar = range(n_steps ** n_params)
        for i in pbar:
            step = [(i // (n_steps ** j)) % n_steps for j in range(n_params)]
            params_i = dict(
                [(param_names[j], 10 ** (param_mins[j] + step[j] * param_steps[j])) for j in range(n_params)])
            self.train(**params_i, **fixed_params)
            cor, gene_cor = self.evaluate()
            if np.mean(cor) > best_perf:
                best_perf = np.mean(cor)
                best_params = params_i

        if stages == 1 or best_params is None:
            if best_params is None:
                best_params = params_i
            train_params = {**best_params, **fixed_params}
            self.train(**train_params)
            return train_params
        else:
            params = dict([(param_names[j], [10 ** (np.log10(best_params[param_names[j]]) - param_steps[j]),
                                             10 ** (np.log10(best_params[param_names[j]]) + param_steps[j])]) for j in
                           range(n_params)])
            return self.grid_search(n_steps, params, fixed_params=fixed_params, verbose=verbose, stages=stages-1, axis=axis)


class Ridge(ModelBase):
    def __init__(self, args, **kwargs):
        super().__init__(args, **kwargs)

    def train(self, fold=0, lam=1):
        if self.svdD[fold] is None:
            self.svdD[fold] = np.linalg.svd(self.D_train[fold], full_matrices=False)
        (ua, sa, va) = self.svdD[fold]
        self.W[fold] = va.T.dot(np.diag(1/(lam*self.n_genes+sa**2)).dot(va)).dot(self.D_train[fold].T.dot(self.Y_train[fold]))


class Stan(ModelBase):
    def __init__(self, args, **kwargs):
        super().__init__(args, **kwargs)

    def train(self, fold=0, lam2=1, lam1=5):
        if self.svdD[fold] is None:
            self.svdD[fold] = np.linalg.svd(self.D_train[fold], full_matrices=False)
        if self.svdK is None:
            uk,sk,vk = np.linalg.svd(self.K, full_matrices=True)
            sk = np.concatenate((sk, [0]*(uk.shape[0]-len(sk))))
            self.svdK = (uk,sk,vk)

        [ud,sd,vd] = self.svdD[fold]
        [uk,sk,vk] = self.svdK
        D = self.D_train[fold]
        Y = self.Y_train[fold]

        # scale = np.divide(1, sd.reshape((-1,1))**2+lw*(1-(lw*sk**2/(lw*sk**2+ls))).reshape(1,-1))
        # W = vd.T.dot(np.multiply(scale, vd.dot(D.T.dot(Y)).dot(uk))).dot(uk.T)
        scale = np.divide(1, sd.reshape((-1,1))**2 + lam1*lam2/(lam1*sk**2+lam2).reshape(1,-1))
        W = vd.T@(np.multiply(scale, vd@D.T@Y@uk))@uk.T
        self.W[fold] = W
