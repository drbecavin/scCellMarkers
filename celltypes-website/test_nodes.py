import load
import yaml
import pandas as pd


# import data for test and data to modify
def import_data_test(selected_compart):
    # Open the yaml file
    with open(
        load.paths(
        )['DATA_PATH'] + f"/yaml_tree/V8/{selected_compart}_tree_formatted_version_V8.yml", 'r') as tree_yaml:
        tree = yaml.load(tree_yaml, Loader=yaml.FullLoader)
    # open the json file to modify
    df_filtered_ref_redcap_sources = pd.read_json(
        load.paths()['DATA_PATH'] + '/results/df_filtered_ref_redcap_sources.json')

    return tree, df_filtered_ref_redcap_sources


def extract_nodes(tree):
    """
    Create a list of nodes in yaml tree
    (keys of nested dictionaries)
    """
    nodes = []
    for k, v in tree.items():
        nodes.append(k)
        if isinstance(v, dict):
            nodes.extend(extract_nodes(v))
        elif isinstance(v, list):
            for item in v:
                if isinstance(item, dict):
                    nodes.extend(extract_nodes(item))
    return nodes


def print_values(search_string, my_dict, all_keys):
    """
    print synonyms of nodes if present
    """

    output = []

    if search_string in all_keys:
        if isinstance(my_dict, list):
            for element in my_dict:
                output += print_values(search_string, element, all_keys)
        elif isinstance(my_dict, dict):
            for key, value in my_dict.items():
                if key == search_string:
                    if isinstance(value, list):
                        for item in value:
                            if isinstance(item, str):
                                print(item)
                                output.append(item)
                    elif isinstance(value, dict):
                        pass
                else:
                    output += print_values(search_string, value, all_keys)
    else:
        pass

    return output


def test_tree_nodes(tree, df_filtered_ref_redcap_sources):
    """
    For a given yaml tree, tests the function 
    """

    # generate the list with the nodes
    all_keys = extract_nodes(tree)

    # print synonyms for all the nodes
    for el in all_keys:
        print('NODE: ', el)
        el_formatted = el.upper()
        el_formatted = el_formatted.strip()
        el_formatted = el_formatted.replace(' ', '__')
        print('NODE formatted is ', el_formatted)
        print('Synonyms are: ')
        synonyms = print_values(el, tree, all_keys)
        synonyms_formatted = []

        for syn in synonyms:
            synonyms_formatted.append(syn.split(';', 1)[0])
        # if the node formatted is not in synonyms
        if el_formatted not in synonyms_formatted:
            print('-----------WARNING!!-----------')
            print('NODE not in synonyms_formatted')
            print('-----------WARNING!!-----------')
            new_row = {'name_dict': 'treearches',
                       'original_key': el,
                       'original_key_formatted': el_formatted,
                       'synonym': None,
                       'synonym_formatted': None,
                       'match_ref': None,
                       'present_in': {}
                       }
            df_filtered_ref_redcap_sources.append(new_row, ignore_index=True)
        print('\n')


selected_compart = 'endothelial'

tree, df_filtered_ref_redcap_sources = import_data_test(selected_compart)

test_tree_nodes(tree)
