import pandas as pd
import load
import search_celltype as search

def search_cell(celltype_search, n):
    df_all = pd.read_json(
        load.paths()['DATA_PATH']+ '/results/df_all.json', 
        orient='records')
    df_filtered_ref_redcap_sources = pd.read_json(
        load.paths()['DATA_PATH']+ '/results/df_filtered_ref_redcap_sources.json', 
        orient='records')
    all_rank_genes = pd.read_json(
        load.paths()['DATA_PATH']+ '/results/all_rank_genes.json', 
        orient='records')

    # Select rows based on the node name
    celltype_search, tmp_df, dict_names = search.search_celltypes_and_dicts(
        df_filtered_ref_redcap_sources, celltype_search)
    # if celltype was found:
    if type(tmp_df) != str:
        dict_genes, rank_genes_selected = search.find_best_genes_for_celltype(celltype_search, df_all,
                                                            tmp_df, dict_names, all_rank_genes, n)
    
        return rank_genes_selected

