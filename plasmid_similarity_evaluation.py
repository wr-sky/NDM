import os


def analyze_similarity():
    limit = 97.5
    blast_result_in = '/home/wbq/4-MCR/2-data_mcr-1/10-plasmid_geography/IncX4/all_97_polished.tsv'
    blast_result_out = '/home/wbq/4-MCR/2-data_mcr-1/10-plasmid_geography/IncX4/all_97_min.tsv'
    results_out = '/home/wbq/4-MCR/2-data_mcr-1/10-plasmid_geography/IncX4/results_limit_{}.tsv'.format(limit)
    basic_file = '/home/wbq/4-MCR/2-data_mcr-1/1-basic_info/basic_info_modified.tsv'

    sp_year = dict()
    sp_prov = dict()
    with open(basic_file, 'r') as fin:
        for line in fin.readlines():
            inf = line.split('\t')
            sp = inf[0]
            prov = inf[1]
            year = inf[2].split('-')[0].split('/')[0]
            if year == '':
                year = '-'
            sp_year[sp] = year
            sp_prov[sp] = prov

    sp_sp_similarity = dict()
    with open(blast_result_in, 'r') as fin:
        for line in fin.readlines():
            inf = line.split('\t')
            for i in range(len(inf)):
                if inf[i][-2] == '_':
                    inf[i] = inf[i].replace('_2', '')
            sp_sp_similarity[inf[0] + '-' + inf[1]] = int(inf[2].replace('\n', ''))

    min_results = dict()
    for key in sp_sp_similarity.keys():
        key_convert = key.split('-')[1] + '-' + key.split('-')[0]
        if not(key in min_results.keys() or key_convert in min_results.keys()):
            count_1 = sp_sp_similarity[key]
            count_2 = sp_sp_similarity[key_convert]
            count = min(count_2, count_1)
            min_results[key] = count

    with open(blast_result_out, 'w') as fout:
        for key in min_results:
            fout.write(key.split('-')[0] + '\t' + key.split('-')[1] + '\t' + str(min_results[key]) + '\n')

    # clustering under limitation

    aa = 1
    sp_sp_list = list()
    for key in min_results:
        if min_results[key] >= limit:
            sp_sp_list.append([key.split('-')[0], key.split('-')[1]])

    all_list = list()
    for key in min_results:
        all_list.append(key.split('-')[0])
        all_list.append(key.split('-')[1])

    cluster_result = list()
    while len(sp_sp_list) > 0:
        base = sp_sp_list[0]
        sp_sp_list.pop(0)
        if len(sp_sp_list) > 0:
            begin = 0
            end = 1
            while begin != end:
                begin = len(sp_sp_list)
                to_remove = list()
                for i in sp_sp_list:
                    if i[0] in base:
                        base.append(i[1])
                        base = list(set(base))
                        to_remove.append(i)
                    elif i[1] in base:
                        base.append(i[0])
                        base = list(set(base))
                        to_remove.append(i)
                for key in to_remove:
                    sp_sp_list.remove(key)
                end = len(sp_sp_list)
        cluster_result.append(base)

    with open(results_out, 'w') as fin:
        for a in cluster_result:
            all_list = list(set(all_list) - set(a))
            fin.write(str(len(a)) + '\n')
            if aa:
                for b in a:
                    fin.write(b + '\t')
                fin.write('\n')
                for b in a:
                    fin.write(sp_prov[b] + '\t')
                fin.write('\n')
                for b in a:
                    fin.write(sp_year[b] + '\t')
                fin.write('\n')
                fin.write('\n')


def cal_gene_sim():
    tsv_path = '/home/wbq/4-MCR/2-data_mcr-1/10-plasmid_geography/IncHI2A/all_147.tsv'
    target_list = ['Salmonella_GCA_015697745.1', 'Salmonella_GCA_023330845.1', 'Salmonella_GCA_023330805.1','Salmonella_GCA_021474285.1']

    # target_list = ['Salmonella_GCA_018531225.1','Salmonella_GCA_018272815.1','Salmonella_GCA_018272585.1','Salmonella_GCA_018272945.1', 'Salmonella_GCA_018272805.1','Salmonella_GCA_018273885.1','Salmonella_GCA_018273915.1','Salmonella_GCA_018272565.1', 'Salmonella_GCA_018272265.1', 'Salmonella_GCA_018273715.1', 'Salmonella_GCA_018273585.1', 'Salmonella_GCA_018272385.1','Salmonella_GCA_023842715.1','Salmonella_GCA_018273565.1','Salmonella_GCA_018273985.1','Salmonella_GCA_018272295.1']
    # target_list = ['Salmonella_GCA_018460045.1', 'Salmonella_GCA_024364965.1', 'Salmonella_GCA_023330765.1', 'Escherichia_GCA_016806045.1', 'Salmonella_GCA_015697805.1', 'Salmonella_GCA_024365025.1']
    # target_list = ['Escherichia_GCA_001894385.1', 'Escherichia_GCA_003858725.1', 'Escherichia_GCA_002941165.1', 'Escherichia_GCA_003859245.1']
    # target_list = ['Salmonella_GCA_030033475.1', 'Escherichia_GCA_003858735.1', 'Escherichia_GCA_018279125.1']
    # target_list = ['Escherichia_GCA_025159115.1', 'Salmonella_GCA_023840185.1', 'Escherichia_GCA_003111785.1', 'Escherichia_GCA_017161345.1', 'Salmonella_GCA_002224405.1', 'Escherichia_GCA_013389695.1', 'Escherichia_GCA_018415995.1', 'Escherichia_GCA_028825485.1', 'Escherichia_GCA_002959275.1', 'Salmonella_GCA_015697825.1', 'Escherichia_GCA_013371745.1', 'Escherichia_GCA_003076395.1']
    # target_list = ['Salmonella_GCA_004124195.2', 'Escherichia_GCA_009909505.1', 'Salmonella_GCA_024364845.1', 'Escherichia_GCA_005233995.1', 'Escherichia_GCA_024499265.1', 'Escherichia_GCA_011492985.1', 'Salmonella_GCA_024364885.1', 'Escherichia_GCA_022559265.1', 'Salmonella_GCA_023330865.1', 'Escherichia_GCA_009659545.1', 'Escherichia_GCA_016446335.1', 'Escherichia_GCA_029026925.1', 'Escherichia_GCA_013201255.1', 'Escherichia_GCA_019375055.2', 'Escherichia_GCA_028622955.1', 'Salmonella_GCA_004010755.1', 'Escherichia_GCA_017868515.1', 'Escherichia_GCA_016889565.1']
    # target_list = ['Escherichia_GCA_019958195.1', 'Escherichia_GCA_019957805.1', 'Escherichia_GCA_019958145.1']
    basic_file = '/home/wbq/4-MCR/2-data_mcr-1/1-basic_info/basic_info_modified.tsv'

    sp_year = dict()
    sp_prov = dict()
    with open(basic_file, 'r') as fin:
        for line in fin.readlines():
            inf = line.split('\t')
            sp = inf[0]
            prov = inf[1]
            year = inf[2].split('-')[0].split('/')[0]
            if year == '':
                year = '-'
            sp_year[sp] = year
            sp_prov[sp] = prov

    sp_sp_score = dict()
    with open(tsv_path, 'r') as fin:
        for line in fin.readlines():
            inf = line.split('\t')
            sp_1 = inf[0]
            sp_2 = inf[1]
            key = sp_1 + '-' + sp_2
            score = float(inf[-2])
            if key in sp_sp_score.keys():
                sp_sp_score[key] = max(sp_sp_score[key], score)
            else:
                sp_sp_score[key] = score

    for sp_1 in target_list:
        tmp = list()
        max_c = 0
        for sp_2 in target_list:
            if sp_2 != sp_1 and sp_year[sp_2] <= sp_year[sp_1]:
                score = sp_sp_score[sp_1+'-'+sp_2]
                if score > max_c:
                    tmp = [sp_2]
                    max_c = score
                elif score == max_c:
                    tmp.append(sp_2)
        print(sp_1 + ' ' + sp_year[sp_1] + ' ' + sp_prov[sp_1] + ': ', end='')
        for a in tmp:
            print(a + ' ' + sp_year[a] + ' ' + sp_prov[a], end=' ')
        print(max_c)


if __name__ == '__main__':
    analyze_similarity()
    cal_gene_sim()
