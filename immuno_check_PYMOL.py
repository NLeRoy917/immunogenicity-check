from pymol import cmd
import requests
import json
import os
from IEDB import *

@cmd.extend
def immuno_check(method='smm',threshold = 1):

    threshold = float(threshold)
    cmd.color('blue','all')

    stored.residues = []
    stored.resvs = []
    stored.bad_resvs = []
    cmd.iterate('all','stored.residues.append(oneletter)')
    cmd.iterate('all','stored.resvs.append(resv)')

    current_resv = 0
    seq = ''

    # Orgaize and clea the sequence
    for resv,aa in zip(stored.resvs,stored.residues):
        if resv == current_resv:
            continue
            # Do nothing and wait for next residue
        else:
            current_resv = resv
            seq += aa

    amino_acids = ['G','P','A','V','L','I','M',
                    'C','F','Y','W','H','K','R',
                    'Q','N','E','D','S','T']
    for aa in seq:
        if aa not in amino_acids:
            #print('Invaid Amino Acid:',aa)
            seq = seq.replace(aa,'')

    flag = 0 # flag for potential immunogenicity detection
    count = 0 # count to skip first line of API response

    # Pack necessary data to send request
    data_pack = {'method':method, 'sequence_text':seq, 'allele':'HLA-A*01:01','length':9}

    #POST request the API and collect the response in Response obj
    try:
        response = requests.post('http://tools-cluster-interface.iedb.org/tools_api/mhci/',data=data_pack)
    except requests.exceptions.RequestException as e:
        print(e)
        sys.stderr.write('Check Network Connection')
        sys.exit(1)

    # Split along new lines
    content = response.text.split('\n')
    first_line = content[0]
    print(first_line)
    if not first_line.startswith('allele'):
        sys.stderr.write(response.text)
        sys.exit(1)

    for line in content[0:-1]:
        #print(line)
        if count == 0:
            count = 1
            continue # skip the first entry

        else:
            data_dump = line.split('\t')

            if float(data_dump[-1]) <= threshold:
                print('--- Potentially immunogenic peptide found! ---')
                stored.bad_resvs.append([data_dump[2],data_dump[3],data_dump[7]])
                print('Start:',data_dump[2])
                print('End:',data_dump[3])
                print('Sequence:',data_dump[5])
                print('IC50 (nMol):',data_dump[6])
                print('% Rank:',data_dump[7])
                print()
                flag = 1

    for coord in stored.bad_resvs:
        #print(coord)
        sele = 'resi' + ' ' + str(coord[0]) + '-' + str(coord[1])
        # cmd.color('red',sele)
        rank = float(coord[2])
        if rank <= 0.1*threshold:
            cmd.color('br9',sele)
        elif rank <= 0.2*threshold:
            cmd.color('br8',sele)
        elif rank <= 0.3*threshold:
            cmd.color('br7',sele)
        elif rank <= 0.4*threshold:
            cmd.color('br6',sele)
        elif rank <= 0.5*threshold:
            cmd.color('br5',sele)
        elif rank <= 0.6*threshold:
            cmd.color('br4',sele)
        elif rank <= 0.7*threshold:
            cmd.color('br3',sele)
        elif rank <= 0.8*threshold:
            cmd.color('br2',sele)
        elif rank <= 0.9*threshold:
            cmd.color('br1',sele)
        else:
            cmd.color('br0',sele)

    #print('Sequence:',seq)
