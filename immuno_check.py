"""
Author: Nathan LeRoy
Date: 6/3/2019

Initial test script to submit requests to the IEDB Analysis Resource database and test
a specific antibody primary sequence for immunogenicity
"""
if __name__ == '__main__':
    import requests
    from Bio import SeqIO
    import json
    import os

    # specify the testing engine/algorithm for the database
    method = 'smm'

    # Iterate through each antibody fasta file in directory
    for fasta_file in os.listdir('./files/fasta'):
        print('Processing file:',fasta_file,'...')

        # Iterate through each sequence/chain in the fasta file
        for record in SeqIO.parse('./files/fasta/'+fasta_file,'fasta'):
            flag = 0 # flag for potential immunogenicity detection
            count = 0 # count to skip first line of API response

            print('Sequence:',record.seq)

            # Pack necessary data to send request
            data_pack = {'method':method, 'sequence_text':str(record.seq), 'allele':'HLA-A*01:01','length':9}

            #POST request the API and collect the response in Response obj
            response = requests.post('http://tools-cluster-interface.iedb.org/tools_api/mhci/',data=data_pack)

            # Split along new lines
            content = response.text.split('\n')

            for line in content[0:-1]:
                #print(line)
                if count == 0:
                    count = 1
                    continue # skip the first entry

                else:
                    data_dump = line.split('\t')

                    if float(data_dump[-1]) <= 1:
                        print('Potentially immunogenic peptide found!')
                        print('Start:',data_dump[2])
                        print('End:',data_dump[3])
                        print('Sequence:',data_dump[5])
                        print('IC50 (nMol):',data_dump[6])
                        print('% Rank:',data_dump[7])
                        print()
                        flag = 1

            if flag == 0:
                print('--- No potentially immunogenic peptides detected ---')
