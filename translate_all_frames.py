#!/bin/bash python3

# Created by Carlos Riera-Ruiz
# Date Jun 27, 2022
# Version 1.0

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Data import CodonTable
import os, shutil, getopt, sys

def usage():
    usage="""
          ####-------------------------------------------------------------------####
          ####-------------------------- Important!! ----------------------------####
          #### Requires Biopython. To install Biopython refer to:                ####
          #### https://biopython.org/wiki/Download                               ####
          #### or                                                                ####
          #### https://anaconda.org/conda-forge/biopython                        ####
          ####                                                                   ####
          #### If you run into problems using translator_all_frames.py or have   ####
          #### any question, contact crieraruiz                                  ####
          ####-------------------------------------------------------------------####
          ####-------------------------------------------------------------------####

          Usage:
          python3 translator_all_frames.py -i <path_to_fasta_file> -m <min_aa_len>

          Options:
              -i (directory)    path to a text file with fasta formated sequences

              -m (integer)      min length of translated ORFs (in amino acids)
                                (do not need to specify, by defauld min_aa_len=10)

              -h                display this help and exit
          
    """
    print(usage)

opts, args = getopt.getopt(sys.argv[1:], 'i:m:h')
# initialize values
min_ORF_len=10
for o,a in opts:
    if o == '-i':
        fasta_file = a     # a text file with fasta formated sequences
    if o == '-m':
        min_ORF_len = int(a)     # a text file with fasta formated sequences
    if o == '-h':
        usage()
        sys.exit(2)     # exits after -h is used

# define codon table
CodonTable = CodonTable.ambiguous_dna_by_id[1]

def make_SeqRecord(seq_id, description, seq):
    '''
    returns a formated sequence in fasta format

    '''
    ## if used stand alone uncoment:
    #from Bio.Seq import Seq
    #from Bio.SeqRecord import SeqRecord

    record = SeqRecord(
        Seq(str(seq)),
        id=seq_id,
        #name=name,
        description=str(description),
    )

    return(record.format("fasta"))

def get_potential_ORFs(right_seq, description, frame):
    '''
    takes as input a sequence of nucleotides and translates it to amino acids
    then splits it using the stop codon "*" as separator
    it also returns the largest sequence between two stop codons
    '''
    # translate sequence
    entry_trans = right_seq.translate('1')
    # get the nt section based on the largest translated ORF
    potential_ORFs = entry_trans.split('*')
    largest_potential_ORF = max(potential_ORFs, key=len)
    frame = frame
    return(right_seq,entry_trans,potential_ORFs, largest_potential_ORF, frame)

def find_met(right_seq,aa_seqs):
    '''
    takes as input a list of potential ORFs and trims their N-terminin until
    it finds the first Methionine
    '''
    ORFs_list=[]
    ORF_list_no_met=[]
    for ORF in aa_seqs:
        if 'M' in ORF:
            index = -1
            for char in ORF:
                index += 1
                if char != 'M':
                    pass
                else:
                    met=ORF[index:]
                    break
            if len(met) > min_ORF_len:
                ORFs_list.append(met)
        else:
            met = ORF
            if len(met) > min_ORF_len:
                ORF_list_no_met.append(met)

    return(ORFs_list, ORF_list_no_met, right_seq)

def get_nt_seq(ORFs,nt_seq,aa_seq):
    # get nt positions from translated ORFs
    nt_seqs_list=[]
    for aa_ORF in ORFs:
        # use first and last five aa of each ORF to spot start and end nt possitions
        ORF_start = aa_ORF[:5]
        ORF_end = aa_ORF[-5:]
        counter = -1
        for i in aa_seq:
            counter += 1
            counter2 = counter + 5
            if str(ORF_start) == str(aa_seq[counter:counter2]):
                ORF1 = int((counter*3))
            if str(ORF_end) == str(aa_seq[counter:counter2]):
                ORF2 = int(counter2*3)

        try:
            # get nt sequences based on the ORF1 and ORF2 possitions
            nt_ORF = nt_seq[ORF1:ORF2]
            nt_seqs_list.append((nt_ORF,ORF1,ORF2,aa_ORF))
        except:
            print('no ORF1 or ORF2')

    return(nt_seqs_list)
 
# remove old files
files_list=['raw_translations.aa','ORFs.aa','ORFs.nt','ORFs_no_met.aa','ORFs_no_met.nt']
def remove_old_files(files_list):
    for file_name in files_list:
        if os.path.isfile(file_name):
            os.remove(file_name)

# remove old files / create files / append sequences
remove_old_files(files_list)

# parse file with fasta file sequences
with open(fasta_file,'r') as fasta:
    handle=SeqIO.parse(fasta,format='fasta')
    # create six coding frames
    for entry in handle:
        description=entry.description
        for frame in 1,2,3,-1,-2,-3:
            #frame=entry.description.split('|')[-1]
            if frame == 1:
                right_seq = entry.seq
            elif frame == 2:
                right_seq = entry.seq[1:]
            elif frame == 3:
                right_seq = entry.seq[2:]
            elif frame == -1:
                right_seq = entry.seq.reverse_complement()
            elif frame == -2:
                right_seq = entry.seq.reverse_complement()[1:]
            elif frame == -3:
                right_seq = entry.seq.reverse_complement()[2:]

            # (entry, entry_trans, potential_ORFs, largest_potential_ORF, frame)
            potential_ORFs=get_potential_ORFs(right_seq, description, frame)

            # go through potential ORFs and keep those just with met, then
            # trim the N-termini of each ORF until reaching the Methioine
            # (ORFs_list, ORF_list_no_met)
            ORFs=find_met(potential_ORFs[0],potential_ORFs[2])

            # find the nucleotide sequence of each found ORF
            ORF_nt=get_nt_seq(ORFs[0],potential_ORFs[0],potential_ORFs[1])
            ORF_nt_no_met=get_nt_seq(ORFs[1],potential_ORFs[0],potential_ORFs[1])

            desc = description + str('_frame_') + str(frame)

            with open('raw_translations.aa', 'a') as files_aa:
                files_aa.write(make_SeqRecord(description,desc,potential_ORFs[1]))

            entry_counter = 0
            for entry_ in ORF_nt:
                entry_counter += 1
                # ORF.aa
                desc2 = desc + str('_ORF_') + str(entry_counter)
                with open('ORFs.aa', 'a') as files_aa:
                    files_aa.write(make_SeqRecord(description,desc2,entry_[3]))

                # ORF.nt
                desc3 = description + str('_frame_') + str(frame) \
                    + str('_ORF_') + str(entry_counter) + str('_from_pos_') \
                    + str(entry_[1]) + str('_to_') + str(entry_[2])
                with open('ORFs.nt', 'a') as files_nt:
                    files_nt.write(make_SeqRecord(description,desc3,entry_[0]))

            entry_counter = 0
            for entry_ in ORF_nt_no_met:
                entry_counter += 1
                # ORF.aa
                desc2 = desc + str('_ORF_') + str(entry_counter)
                with open('ORFs_no_met.aa', 'a') as files_aa:
                    files_aa.write(make_SeqRecord(description,desc2,entry_[3]))

                # ORF.nt
                desc3 = description + str('_frame_') + str(frame) \
                    + str('_ORF_') + str(entry_counter) + str('_from_pos_') \
                    + str(entry_[1]) + str('_to_') + str(entry_[2])
                with open('ORFs_no_met.nt', 'a') as files_nt:
                    files_nt.write(make_SeqRecord(description,desc3,entry_[0]))


