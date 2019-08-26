import requests
from Bio import SeqIO
import json
import os

def checkMHCi(sequence,method='smm',allele='HLA-A*01:01',length=9):
    data_pack = {'method':method, 'sequence_text':sequence, 'allele':allele, 'length':length}

    #POST request the API and collect the response in Response obj
    response = requests.post('http://tools-cluster-interface.iedb.org/tools_api/mhci/',data=data_pack)

    #print(response.text)

    # Split along new lines
    content = response.text.split('\n')

    return content


def checkMHCii(sequence,method='NetMHCIIpan',allele='HLA-DRB1*01:01'):

    # Pack data into dictionary to send with POST request
    data_pack = {'method':method, 'sequence_text':sequence, 'allele':allele}

    #POST request the API and collect the response in Response obj
    response = requests.post('http://tools-cluster-interface.iedb.org/tools_api/mhcii/',data=data_pack)

    #print(response.text)

    # Split along new lines
    content = response.text.split('\n')

    return content


def MHCcomp(MHCi,MHCii,threshold):

    #print('Comparing: {}({}) to {}({})'.format(MHCi[5],MHCi[-1],MHCii[5],MHCii[-1]))

    if ((float(MHCi[-1]) > threshold) or (float(MHCii[-1]) > threshold)):
        # Ignore since threshod of one doesnt match
        return False

    if ((MHCi[5] in MHCii[5]) or (MHCii[5] in MHCi[5])):
        # If the sequences overlap a lot, it is a good candidate for an immunogenic peptide
        print('{}({}) and {}({})'.format(MHCi[5],MHCi[-1],MHCii[5],MHCii[-1]))
        print('Match!')
        return True

    else:
        # Otherwise they dont overlap
        return False


"""##########################################################################################

                                        MAIN SCRIPT

#############################################################################################"""

if __name__ == '__main__':
        # specify the testing engine/algorithm for the database
        threshold = 10

        # Iterate through each antibody fasta file in directory
        for fasta_file in os.listdir('./files/fasta'):
            print('\nProcessing file:',fasta_file,'...')

            # Iterate through each sequence/chain in the fasta file
            for record in SeqIO.parse('./files/fasta/'+fasta_file,'fasta'):

                flag = 0 # flag for potential immunogenicity detection
                count = 0 # count to skip first line of API response
                PIPstore = []

                print('Sequence:',record.seq)

                responseMHCi = checkMHCi(str(record.seq))
                responseMHCii = checkMHCii(str(record.seq))

                for MHCi in responseMHCi[1:-1]:

                    MHCi = MHCi.split('\t')

                    if float(MHCi[-1]) > threshold:
                        #print('Skipping MHCi peptide above threshold')
                        continue

                    for MHCii in responseMHCii[1:-1]:

                        MHCii = MHCii.split('\t')

                        if float(MHCii[-1]) > threshold:
                            #print('Skipping MHCii peptide above threshold')
                            continue

                        if MHCcomp(MHCi,MHCii,threshold):
                            # Store in array
                            PIPstore.append([MHCi,MHCii])
                        else:
                            # Skip and move on
                            continue

                for hit in PIPstore:
                    print('--- Potentially Immunogenic Peptide Found ---')
                    start = hit[0][2]
                    end = hit[1][2]
                    seq = hit[0][5]
                    ic50 = (float(hit[0][6]) + float(hit[1][6]))/2
                    rank_perc = (float(hit[0][7]) + float(hit[1][7]))/2

                    print('Start:',start)
                    print('End:',end)
                    print('Sequence:',seq)
                    print('IC50 (nMol):',ic50)
                    print('% Rank:',rank_perc)
                    print()
