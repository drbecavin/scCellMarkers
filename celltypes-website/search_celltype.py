import pandas as pd

from load import choose_source as source_param


def format_rank_genes_and_sort(rank_genes_df, source, level):
    """
    This function formats a rank_genes_df and sort values
    based on logfold change and p_value
    :param source:
    :param level:
    :param rank_genes_df:
    :return:
    """

    # format data
    rank_genes_df.rename(columns={'group': 'celltype', 'names': 'gene'}, inplace=True)
    rank_genes_df.drop(columns='scores', inplace=True)
    rank_genes_df['source'] = source
    rank_genes_df['level_of_annotation'] = level
    
    # sort based on columns logfoldchanges and pvals_adj
    # we want to keep the same order of the celltypes, but we want to sort
    # columns by not ascending for columns logfoldchanges and pvals_adj
    rank_genes_df_formatted = rank_genes_df.sort_values(['celltype', 'logfoldchanges', 'pvals_adj'],
                                        ascending=[True, False, False])
    rank_genes_df_formatted['idx'] = rank_genes_df.groupby('celltype').cumcount()
    rank_genes_df_formatted = rank_genes_df_formatted.sort_values(['celltype', 'idx'])
    rank_genes_df_formatted = rank_genes_df_formatted.drop(columns=['idx'])
    rank_genes_df_formatted = rank_genes_df_formatted.reset_index(drop=True)
    
    return rank_genes_df_formatted


def format_all_data(rank_genes_dict_GB_Maddi, rank_genes_dict_AJ_Deprez, rank_genes_dict_HLCA,
                    rank_genes_dict_COPD, rank_genes_dict_HTAP):
    """
    This function creates a df with all the rank_genes, each one previously formatted
    and whose values are sorted by using 'logfoldchanges' and 'pvals_adj' columns
    :param rank_genes_dict_GB_Maddi:
    :param rank_genes_dict_AJ_Deprez:
    :param rank_genes_dict_HLCA:
    :param rank_genes_dict_COPD:
    :param rank_genes_dict_HTAP:
    :return:
    """

    # extracting names of levels of annotation for each h5ad file
    name_file_GB_Maddi, groupby_param_GB_Maddi = source_param('GenomeBiology_Maddi')
    name_file_AJ_Deprez, groupby_param_AJ_Deprez = source_param('AJRCCM_Deprez')
    name_file_HLCA, groupby_param_HLCA = source_param('HLCA')
    name_file_COPD, groupby_param_COPD = source_param('COPD')
    name_file_HTAP, groupby_param_HTAP = source_param('HTAP_Yvon')

    # GB_Maddi
    for level in groupby_param_GB_Maddi:
        df_GB_Maddi = format_rank_genes_and_sort(rank_genes_dict_GB_Maddi[level],
                                                 'GenomeBiology_Maddi', level)
    # create final df where we will store all the rank genes information
    # starting by GB_Maddi
    all_rank_genes = df_GB_Maddi.copy()
    
    # AJ_Deprez
    for level in groupby_param_AJ_Deprez:
        df_AJ_Deprez = format_rank_genes_and_sort(rank_genes_dict_AJ_Deprez[level],
                                                  'AJRCCM_Deprez', level)
        # adding rows of AJ_Deprez
        all_rank_genes = pd.concat([all_rank_genes, df_AJ_Deprez])

    # HLCA
    for level in groupby_param_HLCA:
        df_HLCA = format_rank_genes_and_sort(rank_genes_dict_HLCA[level],
                                             'HLCA', level)
        # adding rows of HLCA
        all_rank_genes = pd.concat([all_rank_genes, df_HLCA])
        
    # COPD
    for level in groupby_param_COPD:
        df_COPD = format_rank_genes_and_sort(rank_genes_dict_COPD[level],
                                             'COPD', level)
        # adding rows of COPD
        all_rank_genes = pd.concat([all_rank_genes, df_COPD])

    # HTAP_Yvon
    for level in groupby_param_HTAP:
        df_HTAP = format_rank_genes_and_sort(rank_genes_dict_HTAP[level],
                                             'HTAP_Yvon', level)
        # adding rows of HTAP
        all_rank_genes = pd.concat([all_rank_genes, df_HTAP])
        
    return all_rank_genes


# After all data is formatted, all_rank_genes can be the input of the search_celltype
# search_celltypes_and_dicts
def search_celltypes_and_dicts(df_filtered_ref_redcap_sources, celltype_search):
    """
    This function searches the celltype celltype_search
    in df_filtered_ref_redcap_sources and tells whether it exists or not
    and in which sources.
    :param df_filtered_ref_redcap_sources:
    :param celltype_search:
    :return:
    """

    # names of the dictionaries containing the celltype
    dict_names = []

    # if the celltype is not contained in either original_key or original_key_formatted 
    if (df_filtered_ref_redcap_sources.loc[
        (df_filtered_ref_redcap_sources['original_key'] == celltype_search) | (
                df_filtered_ref_redcap_sources['original_key_formatted'] == celltype_search)].empty):
        print('Sorry, cell type not found')
        tmp_df = 'empty'
    
    # else
    else:
        # creating temporary df after selection of rows
        tmp_df = df_filtered_ref_redcap_sources.loc[
            (df_filtered_ref_redcap_sources['original_key'] == celltype_search) | (
                df_filtered_ref_redcap_sources['original_key_formatted'] == celltype_search)]
        tmp_df.reset_index(inplace=True, drop=True)

        # if column present_in present an empty dict, no synonyms were found
        # hence, only one source is present
        if tmp_df.at[0, 'present_in'] == {}:
            print(celltype_search, ' celltype is present in: ', tmp_df.iloc[0, 0])
        else:
            dict_names.append(tmp_df.iloc[0, 0])
            dict_names += [k for k in tmp_df.at[0, 'present_in'].keys()]
            print(f'celltype {celltype_search} is present in: ', ', '.join(dict_names))
    
    return celltype_search, tmp_df, dict_names


# if celltype was found, it is possible to continue onto the next
# steps: Display top genes in each dict
def print_genes_and_info(all_rank_genes, dict_name, celltype_search, n):
    """
    This function prints information (results of rank_genes) of markers of a given
    celltype_search.
    :param all_rank_genes:
    :param dict_name:
    :param celltype_search:
    :param n:
    :return:
    """
    pd.set_option('display.max_columns', None)
    
    # creating df
    selected_markers_df = all_rank_genes[
        (all_rank_genes['celltype'].str.casefold() == celltype_search.casefold()
        ) & (all_rank_genes['source'].str.casefold() == dict_name.casefold())]

    if (len(selected_markers_df.gene.to_list()) < n) & (len(selected_markers_df.gene.to_list()) > 0):
        print(f'Sorry, less than {n} genes are present...')
        print(f'Best {len(selected_markers_df.gene.to_list())} genes for {celltype_search} in {dict_name} are: \n',
              selected_markers_df.gene.to_list()[:n])
        print('\n \n')
        print('More info about them...')
        selected_markers_df.reset_index(inplace=True)
        print('---------------')
        print('---------------')
        print('\n \n')
    elif (len(selected_markers_df.gene.to_list()) < 1):
        print(f'Sorry, no markers available for {dict_name}')
    else:
        print('\n \n \n')
        print(dict_name.upper(), ": ")
        print(f'Best {n} genes for {celltype_search} in {dict_name} are: \n',
              selected_markers_df.gene.to_list()[:n])
        print('\n \n')
        print('More info about them...')
        selected_markers_df.reset_index(inplace=True)
        print(selected_markers_df.loc[:n])
        print('---------------')
        print('---------------')
        print('\n \n')

    lst_genes = selected_markers_df.gene.to_list()

    return lst_genes, selected_markers_df


def know_original_key_print_info(celltype_search, df_all, all_rank_genes, tmp_df, dict_name, n, dict_genes):
    """
    This function allows to know the original_key of a cell type in a source
    :param celltype_search:
    :param df_all:
    :param all_rank_genes:
    :param tmp_df:
    :param dict_name:
    :param n:
    :param dict_genes:
    :return:
    """
    # the first name
    # This way I know how celltype_search is originally
    # called in the dict
    original_celltype = df_all[
        (df_all['original_key_formatted'] == tmp_df.at[0, "present_in"][dict_name]) & (
            df_all['name_dict'] == dict_name)].iloc[0, 1]
    
    lst_genes, selected_markers_df = print_genes_and_info(all_rank_genes, dict_name, original_celltype, n)
    dict_genes[dict_name] = lst_genes

    return dict_genes, selected_markers_df


def find_best_genes_for_celltype(celltype_search, df_all, tmp_df, dict_names, all_rank_genes, n):
    """

    :param celltype_search:
    :param df_all:
    :param tmp_df:
    :param dict_names:
    :param all_rank_genes:
    :param n:
    :return:
    """
    dict_genes = {}

    #if 'ref_redcap' in dict_names:
    #    dict_names.remove('ref_redcap')

    # Creating a df with all the selected rank genes from all the sources
    rank_genes_selected = pd.DataFrame(columns=all_rank_genes.columns)

    # if the celltype is present only in a dictionary
    # and that dictionary is not ref_redcap or sikkema_discovair or maddi_discovair or deprez_discovair
    if (len(dict_names) == 1) and ((tmp_df.iloc[0, 0] != 'ref_redcap') and (
            tmp_df.iloc[0, 0] != 'sikkema_discovair') and (tmp_df.iloc[0, 0] != 'maddi_discovair') and (
            tmp_df.iloc[0, 0] != 'deprez_discovair')):
        lst_genes, selected_markers_df = print_genes_and_info(all_rank_genes, tmp_df.iloc[0, 0], celltype_search, n)
        rank_genes_selected = pd.concat([rank_genes_selected, selected_markers_df])
        dict_genes[tmp_df.iloc[0, 0]] = lst_genes

    # if the celltype is present only in a dictionary
    # and that dictionary is ref_redcap or sikkema_discovair or maddi_discovair or deprez_discovair
    elif (len(dict_names) == 1) and ((tmp_df.iloc[0, 0] == 'ref_redcap') or (
            tmp_df.iloc[0, 0] == 'sikkema_discovair') or (tmp_df.iloc[0, 0] == 'maddi_discovair') or (
            tmp_df.iloc[0, 0] == 'deprez_discovair')):
        print(f'Sorry... celltype found only in {tmp_df.iloc[0, 0]}... no markers available')

    # if the celltype is present in more than a dictionary
    # and ref_redcap is one of them (in this case it is
    # stored in tmp_df.iloc[0,0]) but sikkema_discovair, maddi_discovair or deprez_discovair are not
    elif (len(dict_names) > 1) and (tmp_df.iloc[0, 0] == 'ref_redcap') and (
            'sikkema_discovair' or 'maddi_discovair' or 'deprez_discovair' not in dict_names):
        dict_names.remove('ref_redcap')
        for dict_name in dict_names:
            dict_genes, selected_markers_df = know_original_key_print_info(celltype_search,
                df_all, all_rank_genes, tmp_df, dict_name, n, dict_genes)
            rank_genes_selected = pd.concat([rank_genes_selected, selected_markers_df])

    # if the celltype is present in more than a dictionary
    # and ref_redcap is not one of them, but
    # sikkema_discovair or maddi_discovair or deprez_discovair is one of them
    elif (len(dict_names) > 1) and (tmp_df.iloc[0, 0] != 'ref_redcap') and (
            'sikkema_discovair' or 'maddi_discovair' or 'deprez_discovair' in dict_names):
        # if the first dict is sikkema_discovair or maddi_discovair or deprez_discovair
        if (tmp_df.iloc[0,0] == 'sikkema_discovair') or (
            tmp_df.iloc[0,0] == 'maddi_discovair') or (
            tmp_df.iloc[0,0] == 'deprez_discovair'):
            dict_names.remove(tmp_df.iloc[0,0])
        else: 
            lst_genes, selected_markers_df = print_genes_and_info(all_rank_genes, tmp_df.iloc[0, 0], celltype_search, n)
            rank_genes_selected = pd.concat([rank_genes_selected, selected_markers_df])
            dict_names.remove(tmp_df.iloc[0,0])
        for dict_name in dict_names:
            if (dict_name != 'sikkema_discovair') and (dict_name != 'maddi_discovair') and (
                    dict_name != 'deprez_discovair'):
                dict_genes, selected_markers_df = know_original_key_print_info(celltype_search, 
                                                                               df_all, all_rank_genes,
                                                                               tmp_df, dict_name, 
                                                                               n, dict_genes)
                rank_genes_selected = pd.concat([rank_genes_selected, selected_markers_df])
            else:
                continue

    # if the celltype is present in more than a dictionary and
    # both ref_redcap and one among sikkema_discovair or maddi_discovair or deprez_discovair is one of them
    elif (len(dict_names) > 1) and (tmp_df.iloc[0, 0] == 'ref_redcap') and (
            'sikkema_discovair' or 'maddi_discovair' or 'deprez_discovair' in dict_names):
        dict_names.remove(tmp_df.iloc[0,0])
        for dict_name in dict_names:
            if (dict_name != 'sikkema_discovair') and (dict_name != 'maddi_discovair') and (
                    dict_name != 'deprez_discovair'):
                dict_genes, selected_markers_df = know_original_key_print_info(celltype_search,
                                                                               df_all, all_rank_genes, 
                                                                               tmp_df, dict_name, n, dict_genes)
                rank_genes_selected = pd.concat([rank_genes_selected, selected_markers_df])
            else:
                continue

    # if the celltype is present in more than a dictionary
    # and ref_redcap is not one of them and
    # sikkema_discovair or maddi_discovair or deprez_discovair as well are not
    elif (len(dict_names) > 1) and (tmp_df.iloc[0, 0] != 'ref_redcap') and (
            'sikkema_discovair' or 'maddi_discovair' or 'deprez_discovair' not in dict_names):
        lst_genes, selected_markers_df = print_genes_and_info(all_rank_genes, tmp_df.iloc[0, 0], celltype_search, n)
        rank_genes_selected = pd.concat([rank_genes_selected, selected_markers_df])
        dict_names.remove(tmp_df.iloc[0,0])
        dict_genes[tmp_df.iloc[0, 0]] = lst_genes
        for dict_name in dict_names:
            dict_genes, selected_markers_df = know_original_key_print_info(
                df_all, all_rank_genes, tmp_df, dict_name, n, dict_genes)
            rank_genes_selected = pd.concat([rank_genes_selected, selected_markers_df])

    return dict_genes, rank_genes_selected
