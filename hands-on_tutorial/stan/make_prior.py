import os
import pandas as pd
import scanpy as sc
from anndata import AnnData
from pandas import DataFrame
from typing import Optional, Union, Literal, List


def make_dir(dir):
    if not os.path.exists(dir):
       # Create a new directory because it does not exist
       os.makedirs(dir)


def load_humantfs(source_dir="data/gene_tf/"):
    fname = source_dir+"humantfs.csv"
    if os.path.isfile(fname):
        humantfs = pd.read_csv(fname)
    else:
        print('Downloading humantfs database...')
        humantfs = pd.read_csv("http://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.csv")
        make_dir(source_dir)
        humantfs.to_csv(fname)
    return humantfs.query("`Is TF?`=='Yes' and `TF assessment` == 'Known motif'")['HGNC symbol'].tolist()


def load_hTFtarget(source_dir="data/gene_tf/"):
    fname = source_dir+"hTFtarget.csv"
    if os.path.isfile(fname):
        htftarget = pd.read_csv(fname, index_col=0)
    else:
        print('Downloading hTFtarget database...')
        hTFtarget_url = "http://bioinfo.life.hust.edu.cn/static/hTFtarget/file_download/tf-target-infomation.txt"
        htftarget = pd.read_csv(hTFtarget_url,sep='\t')
        htftarget.columns = ['TF', 'gene', 'tissue']
        make_dir(source_dir)
        htftarget.to_csv(fname)

    htftarget = htftarget.drop(columns="tissue")
    htftarget.columns = ["TF", "gene"]
    return htftarget


def add_gene_tf_matrix(
        data: AnnData,
        min_cells_proportion: float = 0.2,
        min_tfs_per_gene: int = 5,
        min_genes_per_tf: int = 10,
        gene_tf_source: Union[Literal["dorothea", "hTFtarget"], DataFrame]="hTFtarget",
        tf_list: Union["humantfs", List, None]="humantfs",
        source_dir: str = "data/gene_tf/"
        ) -> Union[AnnData, None]:

    adata = data.copy()
    tf_universe = adata.var_names.copy()
    sc.pp.filter_genes(adata,
                       min_cells=min_cells_proportion*adata.n_obs)

    gene_tf = make_gene_tf_matrix(source=gene_tf_source,
                                tf_list="humantfs",
                                gene_universe=adata.var_names,
                                tf_universe=tf_universe,
                                min_tfs_per_gene=min_tfs_per_gene,
                                min_genes_per_tf=min_genes_per_tf,
                                source_dir=source_dir)

    adata = adata[:, gene_tf.index]
    adata.varm['gene_tf'] = gene_tf
    adata.uns['tf_names'] = gene_tf.columns.to_list()
    return adata
    

def make_gene_tf_matrix(source="hTFtarget",
                        tf_list="humantfs",
                        gene_universe=None,
                        tf_universe=None,
                        min_tfs_per_gene=10,
                        min_genes_per_tf=10,
                        source_dir: str = "data/gene_tf/"):

    if type(source)== pd.DataFrame:
        gene_tf=source
    elif source=="hTFtarget":
        gene_tf=load_hTFtarget(source_dir=source_dir)

    if tf_list=="humantfs":
        tf_list=load_humantfs(source_dir=source_dir)
    if tf_list is not None:
        gene_tf=gene_tf.query("TF in @tf_list")

    if tf_universe is not None:
        gene_tf=gene_tf.query("TF in @tf_universe") #just say >0
    if gene_universe is not None:
        gene_tf=gene_tf.query("gene in @gene_universe")

    gene_tf['values'] = 1
    gene_tf = gene_tf.pivot_table(index='gene', columns='TF', aggfunc='mean', values='values',fill_value=0)

    # linan's version
    gene_sum = gene_tf.T.sum()
    gene_tf.drop(gene_tf.index[gene_sum<min_tfs_per_gene], axis=0,inplace=True)
    tf_sum = gene_tf.sum()
    gene_tf.drop(gene_tf.columns[tf_sum<min_genes_per_tf], axis=1,inplace=True)
    return gene_tf