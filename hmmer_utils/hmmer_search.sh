#!/bin/bash
module load muscle
module load hmmer

# user input
evalue=1e-2
query_name=enzyme
db_path=seqdb.fasta

# HMMER outputs
query_path=${query_name}.fasta
msa_path=outputs/${query_name}_align.fasta
sto_path=outputs/${query_name}_align.sto
hmm_path=outputs/${query_name}.hmm
hmmoutput_path=outputs/${query_name}_hmmeroutput.txt
out_path=outputs/${query_name}_out.id
fasta_path=outputs/${query_name}_${evalue}_hmmsearch.fasta

mkdir outputs

# msa on input sequence fasta
muscle -in ${query_path} -out ${msa_path} -quiet

# convert fasta to stockholm
python ./scripts/fasta2sto.py ${msa_path} ${sto_path}

# build hmm profile for input proteins
hmmbuild ${hmm_path} ${sto_path}

# search in database for homologous proteins
hmmsearch --tblout ${hmmoutput_path} -E ${evalue} --seed 17 ${hmm_path} ${db_path}

# extract ids of homologous proteins
python ./scripts/extractid.py ${hmmoutput_path} ${out_path}

# output fasta of homologous proteins found by hmmer
python ./scripts/extractfasta.py ${db_path} ${out_path} ${fasta_path}

