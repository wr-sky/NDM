import os


def main():
    file_gb = "path_to_gbff_files"
    file_record = "path_to_tsv_files_NDM_downstream/upstream"
    file_filter = "path_to_formatted_tsv_files_NDM_downstream/upstream"

    count = 0
    with open(file_record, 'w') as fout:
        with open(file_gb, 'r') as fin:
            name = ''
            gen_list = list()
            for line in fin.readlines():
                if 'LOCUS' in line:
                    name += line.split('       ')[1]
                if 'DEFINITION' in line:
                    name += '_' + line.split('  ')[1].split(' ')[0]
                # extract gene name
                if '/product=' in line:
                    if len(line.split('\"')) < 2:
                        print(name+'!!!!')
                        print(line)
                    gen_list.append(line.split('\"')[1])
                if '/' == line[0]:
                    if len(gen_list) < 10:
                        count += 1
                    else:
                        fout.write('>' + name + '\t')
                        position = -1
                        flag = False

                        # to detect the direction of gene sequence (forward or reverse)
                        for gen in gen_list:
                            position += 1
                            if 'NDM' in gen or 'Delhi' in gen or 'Ndm' in gen or 'Metallo-beta-lactamase' in gen or 'metallo-beta-lactamase' in gen or ('beta-lactamase' in gen and 'B1' in gen):
                                flag = True
                                break
                        constant = len(gen_list)

                        if not flag:
                            print(name.split('_')[0])
                        else:
                            # upstream_2, downstream_7
                            if 'ble' in gen_list[(position+1) % constant] or 'MBL' in gen_list[(position+1) % constant] or 'Ble' in gen_list[(position+1) % constant] or 'isomerase' in gen_list[(position+2) % constant] or 'trpF' in gen_list[(position+1) % constant]:
                                # Forward sequence
                                for i in range(-2, 10):
                                    fout.write(gen_list[(position+i) % constant].replace('\n', ''))
                                    fout.write('\t')
                                fout.write('\n')
                            else:
                                # reverse sequence
                                for i in range(-2, 10):
                                     fout.write(gen_list[(position-i)%constant].replace('\n', ''))
                                     fout.write('\t')
                                fout.write('\n')
                    name = ''
                    gen_list.clear()
            print('count = ' + str(count))

    # format gene names
    with open(file_record, 'r') as fin:
        with open(file_filter, 'w') as fout:
            for line in fin.readlines():
                if '>' not in line:
                    for gen in line.split('\t\t'):
                        if 'NDM' in gen or 'beta' in gen or 'Delhi' in gen or 'Ndm' in gen:
                            fout.write('NDM-5\t\t')
                            continue
                        if 'bleomycin' in gen or 'bleMBL' in gen or 'MBL' in gen or 'Bleomycin' in gen or 'BleMLB' in gen:
                            fout.write('ble\t\t')
                            continue
                        if 'hypothetical' in gen:
                            fout.write('hypothetical\t\t')
                            continue
                        if 'isomerase' in gen or 'TrpF' in gen:
                            fout.write('trpF\t\t')
                            continue
                        if 'TAT' in gen or 'Dsbc' in gen or 'DsbD' in gen or 'DsbC' in gen:
                            fout.write('tat\t\t')
                            continue
                        if 'tolerance' in gen:
                            fout.write('dct\t\t')
                            continue
                        if 'InsH' in gen:
                            fout.write('InsH\t\t')
                            continue
                        if 'obile' in gen:
                            fout.write('mobile\t\t')
                            continue
                        if 'family' in gen:
                            fout.write(gen.split()[0] + '\t\t')
                            continue
                        fout.write(gen.replace('\n', '') + '\t\t')
                    fout.write('\n')
                else:
                    fout.write(line)


if __name__ == '__main__':
    main()