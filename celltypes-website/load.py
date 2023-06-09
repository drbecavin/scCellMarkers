import os
import re

import pandas as pd
import graphviz
import scanpy as sc


# import_data
# This module contains the functions to upload different
# kind of input data.


def paths():
    """
    This function creates a dict. The values are addresses.
    The key specifies which is the address or file.
    The function returns the dict.

    :return:
    """
    
    current_directory = os.getcwd()
    working_directory = os.path.dirname(current_directory)
    return ({'WD_PATH': working_directory,
             'CD_PATH': current_directory,
             'DATA_PATH': os.path.join(current_directory, 'data')})


# Ref_redcap
def upload_ref_redcap():
    """
    Uploading reference of cell types redcap.
    :return:
    """
    # Importing the reference redcap
    redcap_file = pd.ExcelFile(paths()[
                                   'DATA_PATH'] + '/ref/Marqueurs_RedCAP.xlsx')
    # df1 = pd.read_excel(redcap_file, 'Feuil1')
    redcap_raw = pd.read_excel(redcap_file, "Feuil2")

    # filtering human cell types
    redcap_df = redcap_raw[redcap_raw['Unnamed: 0'].str.contains('Human')]

    # Are there any duplicates?
    # len(redcap_df['Unnamed: 0'].unique()) == redcap_df.shape[0]

    # keeping only columns with the cell types and their synonyms
    redcap_df.drop(redcap_df.columns[[0, 1, 2, 4, 5]],
                   axis=1,
                   inplace=True)
    redcap_df.reset_index(
        drop=True, inplace=True)
    # replacing NaN with None
    redcap_df.where(
        pd.notnull(redcap_df), None,
        inplace=True)

    # creating dictionary where all the cell types
    # and the synonyms are stored
    redcap = {}
    # iterating over each cell type
    for i in range(redcap_df.shape[0]):
        redcap[redcap_df.iloc[i, 0]] = []
        # iterating over all the cells at
        # that row (except for the cell type cell)
        for j in range(1, redcap_df.shape[1]):
            if redcap_df.iloc[i, j] is not None:
                redcap[redcap_df.iloc[i, 0]].append(
                    redcap_df.iloc[i, j])
    return redcap


# Trees in dot format
def upload_dot_tree(cell_group='', version=''):
    """
    This function uploads dot files
    containing trees and creates a Graph object.
    :param cell_group:
    :param version:
    :return:
    """
    # Load the .dot file and create a Graph object

    # if no save name and no version are specified
    if (cell_group == "") & (version == ""):
        dot_path = paths()['DATA_PATH'] + f'/graphviz_tree/trees/Digraph.gv.dot'

    # if a cell_group is specified, but no version is
    elif (cell_group != "") & (version == ''):
        dot_path = paths()['DATA_PATH'] + f'/graphviz_tree/trees/{cell_group}/Digraph.gv.dot'

    # if both cell_group and version are specified
    else:
        dot_path = paths()['DATA_PATH'] + f'/graphviz_tree/trees_V{version}/{cell_group}/Digraph.gv.dot'

    graph = graphviz.Source.from_file(dot_path)
    return graph


def choose_source(source):
    """
    This function allows to choose a source and automatically store parameters.
    It will be used within next functions not to repeat it.
    :param source:
    :return:
    """

    if source == "GenomeBiology_Maddi":
        name_file = paths()['DATA_PATH'] + "/h5ad/madissoon19_lung.processed.h5ad"
        groupby_param = ["CellType"]
    elif source == "AJRCCM_Deprez":
        name_file = paths()['DATA_PATH'] + "/h5ad/deprez19_restricted.processed.h5ad"
        groupby_param = ["CellType"]
    # HLCA, sikkema
    elif source == "HLCA":
        name_file = paths()['DATA_PATH'] + "/h5ad/HLCA_v1.h5ad"
        groupby_param = ["ann_level_1", "ann_level_2", "ann_level_3", "ann_level_4", "ann_level_5"]
    # COPD, collin
    elif source == "COPD":
        name_file = paths()['DATA_PATH'] + "/h5ad/integrated_V9.h5ad"
        groupby_param = ["celltype_lv2_V4", "celltype_lv1_V4", "celltype_lv0_V4"]
    elif source == 'HTAP_Yvon':
        name_file = paths()['DATA_PATH'] + "/h5ad/htap.h5ad"
        groupby_param = ['population', 'lineage', 'celltype']

    return name_file, groupby_param


# rank genes function, it will be an option to use it or not
# in upload_h5ad_adata
def rank_genes(adata, level, source):
    """
    This function run the following two functions: sc.tl.rank_genes_groups and sc.get.rank_genes_groups_df
    on a given adata object, with given parameters: groupby_param and source
    :param level:
    :param adata:
    :param source:
    :return:
    """
    # Calculating rank of genes
    sc.tl.rank_genes_groups(
        adata,  # data
        groupby=level,  # how to group cells (lines) can change from source to source
        method='wilcoxon',  # DA method
        key_added=f'DE_CellType_{level}',  # where are results stored in adata object
        use_raw=False,  # whether is raw data or not.
        pts=True
    )

    # save results to df and to csv
    sc.get.rank_genes_groups_df(adata,
                                None,  # group to return results from.
                                key=f'DE_CellType_{level}'  # key is the argument to choose where to get data
                                ).to_csv(paths()['DATA_PATH'] + f"/rank_genes/rank_genes_{source}_{level}.csv")


# table with rank_genes from Scanpy
def upload_scanpy_rank_genes(source, level):
    """
    This function reads scanpy
    :param level:
    :param source:
    :return:
    """
    return pd.read_csv(paths()['DATA_PATH'] + f"/rank_genes/rank_genes_{source}_{level}.csv", index_col=0)


# .h5ad files
def upload_h5ad_adata(source, rank_genes_bool=False):
    """
    This function reads a h5ad file, and if rank_genes,
    it also calculates rank genes groups and saves to a csv file.
    :param source:
    :param rank_genes_bool:
    :return:
    """

    # choosing source and getting its params
    name_file, groupby_param = choose_source(source)

    # reading h5ad file to create adata object
    adata = sc.read_h5ad(name_file)

    # if user wants to calculate rank_genes
    if rank_genes_bool:
        # creating output folder if it does not exist
        if not os.path.exists(paths()['DATA_PATH'] + '/rank_genes'):
            os.mkdir(paths()['DATA_PATH'] + '/rank_genes')
        # computing rank_genes for all levels of annotation
        for level in groupby_param:
            rank_genes(adata, level, source)

    # Whether it was calculated using this function or not
    # with every .h5ad file it will come also a rank genes df
    rank_genes_dict = {}
    lst_celltypes_all_levels = []
    dict_celltypes = {}
    for level in groupby_param:
        rank_genes_dict[level] = upload_scanpy_rank_genes(source, level)
        # adding celltypes of a source to a list
        lst_celltypes = adata.obs[level].unique().to_list()
        for el in lst_celltypes:
            if el not in lst_celltypes_all_levels:
                lst_celltypes_all_levels.append(el)
    # storing information of celltypes in a dictionary.
    # We need to store information in this format to later
    # use it in the preprocessing part that takes the cell types
    # input as keys of a dictionary.
    for el in lst_celltypes_all_levels:
        dict_celltypes[el] = None

    return adata, rank_genes_dict, dict_celltypes


# Specific cell types and markers from Antoine's work
def upload_collin_other(source):
    """
    This function uploads data from different sources
    of Antoine Collin's work combined with different sources
    :param source:
    :return:
    """
    if not os.path.exists(paths()['DATA_PATH'] + '/collin_other'):
        os.mkdir(paths()['DATA_PATH'] + '/collin_other')

    # sikkema
    if source == 'sikkema_discovair':
        # last version of markers given by Pascal Barbry 20/02/2023
        # upload file
        markers_sikkema_allinfo = pd.read_excel(
            paths()['DATA_PATH'] + '/collin_other/146445_2_supp_1184000_rpwm87.xlsx',
            index_col=0, header=1)

        # filter columns with markers
        markers_sikkema_df = markers_sikkema_allinfo.filter(
            regex='marker$', axis=1)

        # create dictionary to return
        markers = {}
        # for each cell type
        for col in markers_sikkema_df:
            # creating a list of marker genes for cell type col
            lst_mark = markers_sikkema_df.loc[:, col]
            # for each element of the list
            col = re.sub('_marker', "", col)
            # removing Nan
            markers[col] = [x for x in lst_mark if str(x) != 'nan']

    # maddi
    elif source == 'maddi_discovair':
        # upload file
        markers_madi = pd.read_csv(paths()['DATA_PATH'] + '/collin_other/Markers_Maddison.csv', sep=';', index_col=0)
        # formatting data
        markers_madi = markers_madi.fillna(0)
        markers_madi = markers_madi.astype('int').astype('bool')
        markers = {}
        for celltype in markers_madi.columns:
            markers[celltype] = list(markers_madi.loc[:, celltype][markers_madi.loc[:, celltype]].index)
            markers[celltype] = list(pd.Series(markers[celltype]).drop_duplicates())

    # deprez
    elif source == 'deprez_discovair':
        markers_Deprez_df = pd.read_csv(paths()['DATA_PATH'] + '/collin_other/MarkerCells_Deprez.csv', sep=';')
        markers = dict.fromkeys(markers_Deprez_df.columns)
        for col in markers_Deprez_df.columns:
            markers[col] = list([gene for gene in markers_Deprez_df[col].dropna()])

    return markers


def import_all(rank_genes_bool=False, version_trees='8'):
    """
    Final function of module load.
    This function automatically uploads all data.
    If rank_genes_bool: it also calculates rank genes groups and saves to a csv file.
    :param version_trees:
    :param rank_genes_bool:
    :return:
    """
    # Collin other celltypes data
    # sikkema_discovair
    sikkema_discovair = upload_collin_other('sikkema_discovair')
    # maddi_discovair
    maddi_discovair = upload_collin_other('maddi_discovair')
    # deprez_discovair
    deprez_discovair = upload_collin_other('deprez_discovair')

    # Redcap
    ref_redcap = upload_ref_redcap()

    # Trees
    # endothelial
    graph_endothelial = upload_dot_tree(cell_group='endothelial', version=version_trees)
    # epithelial
    graph_epithelial = upload_dot_tree(cell_group='epithelial', version=version_trees)
    # stromal
    graph_stromal = upload_dot_tree(cell_group='stromal', version=version_trees)
    # immune
    graph_immune = upload_dot_tree(cell_group='immune', version=version_trees)

    # h5ad data
    # GenomeBiology_Maddi => GB_Maddi
    adata_GB_Maddi, rank_genes_dict_GB_Maddi, dict_celltypes_GB_Maddi = upload_h5ad_adata(
        source="GenomeBiology_Maddi", rank_genes_bool=rank_genes_bool)
    # AJRCCM_Deprez => AJ_Deprez
    adata_AJ_Deprez, rank_genes_dict_AJ_Deprez, dict_celltypes_AJ_Deprez = upload_h5ad_adata(
        source="AJRCCM_Deprez", rank_genes_bool=rank_genes_bool)
    # HLCA
    adata_HLCA, rank_genes_dict_HLCA, dict_celltypes_HLCA = upload_h5ad_adata(
        "HLCA", rank_genes_bool=rank_genes_bool)
    # COPD
    adata_COPD, rank_genes_dict_COPD, dict_celltypes_COPD = upload_h5ad_adata(
        "COPD", rank_genes_bool=rank_genes_bool)
    # HTAP_Yvon => HTAP
    adata_HTAP, rank_genes_dict_HTAP, dict_celltypes_HTAP = upload_h5ad_adata(
        "HTAP_Yvon", rank_genes_bool=rank_genes_bool)

    return adata_GB_Maddi, rank_genes_dict_GB_Maddi, dict_celltypes_GB_Maddi, \
        adata_AJ_Deprez, rank_genes_dict_AJ_Deprez, dict_celltypes_AJ_Deprez, \
        adata_HLCA, rank_genes_dict_HLCA, dict_celltypes_HLCA, \
        adata_COPD, rank_genes_dict_COPD, dict_celltypes_COPD, \
        adata_HTAP, rank_genes_dict_HTAP, dict_celltypes_HTAP, \
        sikkema_discovair, maddi_discovair, deprez_discovair, \
        ref_redcap, \
        graph_endothelial, graph_epithelial, graph_stromal, graph_immune