import networkx as nx
import os


def cal_score(list_1, list_2):
    score = 0
    # same position and same gene (exclude // and --) +1 +2 +3 +...
    iter_count = 1
    for i in range(len(list_1)):
        if list_1[i] == list_2[i] and (list_1[i] != '--' or list_1[i] != '//'):
            score += iter_count
            iter_count += 1
        else:
            break

    # number of the same gene (exclude // and --) + number
    if len(list_1) > 0:
        set_1 = set(list_1)
        set_1.discard('--')
        set_1.discard('//')
        set_2 = set(list_2)
        set_2.discard('--')
        set_2.discard('//')
        for a in list(set_1.intersection(set_2)):
            score += min(list_1.count(a), list_2.count(a))
    return score


def iter_search(node, temp_dict, temp_list):
    for key in temp_dict.keys():
        if node in key:
            keys = key.split('-')
            keys.remove(node)
            node_2 = keys[0]
            if node_2 not in temp_list:
                temp_list.append(node_2)
                iter_search(node_2, temp_dict, temp_list)


def data_process():
    base_path = 'base_path_to_NDM_down/upstream_genes_tsv_file'
    path = os.path.join(base_path, 'path_to_NDM_down/upstream_genes_tsv_file')

    # read tsv and preprocess data
    inf_dict = dict()
    with open(path, 'r') as fin:
        for line in fin.readlines():
            inf_list = line.replace('\n', '').split('\t')
            inf_list_pro = list()
            for elm in inf_list:
                elm = elm.split(' (')[0]
                # transfer all IS,Tn to transposase
                if 'IS' in elm or 'Tn' in elm:
                    elm = 'transposase'
                # transfer all *S ribosomal RNA to ribosomal
                if 'ribosomal' in elm:
                    elm = 'ribosomal'
                inf_list_pro.append(elm)
            inf_dict[inf_list_pro[0]] = inf_list_pro[1:]

    # calculate score of each pair of sequences
    similarity_dict = dict()
    key_list = sorted(inf_dict.keys())
    for key_1 in key_list:
        for key_2 in key_list:
            if key_list.index(key_1) < key_list.index(key_2):
                list_1 = inf_dict[key_1]
                list_2 = inf_dict[key_2]
                s = cal_score(list_1, list_2)
                similarity_dict[key_1 + '-' + key_2] = s

    # clustering according to the score/fully connected network and record the tree structure into multiple tsv files.
    limit = 35
    node_list = key_list
    temp_dict = similarity_dict

    key_to_del_list = list()
    for key in temp_dict.keys():
        if temp_dict[key] < limit:
            key_to_del_list.append(key)
    for key_to_del in key_to_del_list:
        del temp_dict[key_to_del]
    print(len(temp_dict.keys()))

    list_of_node_list = list()
    while len(node_list) > 0:
        node = node_list[0]
        temp_list = list()
        temp_list.append(node)
        iter_search(node, temp_dict, temp_list)
        list_of_node_list.append(temp_list)
        for temp_node in temp_list:
            node_list.remove(temp_node)
    for list_of_node in list_of_node_list:
        print(list_of_node)

    path_2 = os.path.join(base_path, 'path_to_the_trunk_tsv'.format(str(limit)))
    path_3 = os.path.join(base_path, 'path_to_the_branch_tsv'.format(str(limit)))
    count = 0
    with open(path_2, 'w') as fin:
        for list_of_node in list_of_node_list:
            if len(list_of_node) < 10:
                for name in list_of_node:
                    fin.write(name)
                    for value in inf_dict[name]:
                        fin.write('\t')
                        fin.write(value)
                    fin.write('\n')
            else:
                with open(path_3.replace('.tsv', '_'+str(count)+'.tsv'), 'w') as fin_2:
                    for name in list_of_node:
                        fin_2.write(name)
                        for value in inf_dict[name]:
                            fin_2.write('\t')
                            fin_2.write(value)
                        fin_2.write('\n')
                    count += 1


if __name__ == '__main__':
    data_process()