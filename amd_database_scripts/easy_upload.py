#!/usr/bin/env python
'''
Matt Lueder

This script uploads a phage genome, and uses the genome as input into glimmer3 and trnascan-se to automatically add
cds and trna features. It wraps the upload_phage.py and create_features_for_phage.py scripts. Annotations for tRNAs
and terminal repeats are created automatically.

Requires the tools long-orfs, extract, build-icm, glimmer3, and trnascan-se to be in PATH.
These are available as conda packages:
> conda install glimmer trnascan-se

To see usage information use the option '-h'.
'''

import subprocess
import sys
from tempfile import TemporaryDirectory
import os, django

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "LIMS.settings")
django.setup()


def run_glimmer(fasta_path, phage_name, output_destination, minimum_gene_len=90):
    '''
    Run glimmer and return path to predict file.
    '''
    print('*** RUNNING GLIMMER ***')

    # Run long-orfs
    long_orf_file = open(os.path.join(output_destination, '%s.glimmer.longorfs' % phage_name), 'w')
    subprocess.run(['long-orfs', '-n', '-t', '1.1', fasta_path, '-'], stdout=long_orf_file, check=True)

    # Run extract
    extract_file = open(os.path.join(output_destination, '%s.glimmer.train' % phage_name), 'w+')
    subprocess.run(['extract', '-t', fasta_path, long_orf_file.name], stdout=extract_file, check=True)
    extract_file.seek(0)

    # Run build icm
    icm_file_path = os.path.join(output_destination, '%s.glimmer.icm' % phage_name)
    subprocess.run(['build-icm', '-r', icm_file_path], check=True, stdin=extract_file)

    # Run glimmer3
    subprocess.run(['glimmer3', '-o50', '-g{}'.format(minimum_gene_len), '-t30', fasta_path, icm_file_path,
                    os.path.join(output_destination, '%s.glimmer' % phage_name)], check=True)

    long_orf_file.close()
    extract_file.close()
    # Unless I am mistaken, it may be worth while to add to this the second round of glimmer;
    # I think this, as written, is intended to provide the training hmm for a more refined
    # search by glimmer, as per Arthur Delcher's release notes:
    # https://ccb.jhu.edu/software/glimmer/glim302notes.pdf
    #
    # Thus I think one would add a variant of the end of 'g3-iterated.csh', e.g.:
    # tail +2 run3.run1.predict > run3.coords
    # upstream-coords.awk 25 0 run3.coords | extract genom.seq - > run3.upstream
    # elph run3.upstream LEN=6 | get-motif-counts.awk > run3.motif
    # set startuse = $(start-codon-distrib -3 genom.seq run3.coords)
    # glimmer3 -o50 -g110 -t30 -b run3.motif -P ${startuse} genom.seq run3.icm run3
    #
    # My own script wraps the two iterations in a set of invocations of all the pieces
    # which when combined, look something like this:
    #
    # mkdir -p outputs/glimmer
    # long-orfs -n -t 1.15 outputs/shovill_r1_trimmed-corrected/contigs.fa outputs/glimmer/first_run_longorfs.txt
    # extract outputs/shovill_r1_trimmed-corrected/contigs.fa outputs/glimmer/first_run_longorfs.txt \
    #    1>outputs/glimmer/first_run_training.txt 2>outputs/glimmer/first_run_training.err
    # build-icm -r outputs/glimmer/first_run.icm < outputs/glimmer/first_run_training.txt \
    #    2>outputs/glimmer/first_run_icm.out 1>&2
    # glimmer3 -o50 -g110 -t30 \
    #    outputs/shovill_r1_trimmed-corrected/contigs.fa \
    #    outputs/glimmer/first_run.icm \
    #    outputs/glimmer/first_run.out
    # ##  Use this first run output to go again
    # extract -t outputs/shovill_r1_trimmed-corrected/contigs.fa \
    #    outputs/glimmer/first_run.out.predict train.coords \
    #    1>outputs/glimmer/second_run_training.txt \
    #    2>outputs/glimmer/second_run_training.err
    # build-icm -r outputs/glimmer/second_run.icm < outputs/glimmer/second_run_training.txt
    # upstream-coords.awk 25 0 outputs/glimmer/first_run.out.predict | extract outputs/shovill_r1_trimmed-corrected/contigs.fa - > \
    #    outputs/glimmer/second_run_upstream.txt
    # elph outputs/glimmer/second_run_upstream.txt LEN=6 2>outputs/glimmer/elph.err | \
    #    get-motif-counts.awk 2>outputs/glimmer/second_run_motif.txt 1>&2
    # startuse=$(start-codon-distrib -3 outputs/shovill_r1_trimmed-corrected/contigs.fa outputs/glimmer/first_run.out.predict)
    # glimmer3 -o50 -g110 -t30 -b outputs/glimmer/second_run_motif.txt -P ${startuse} \
    #    outputs/shovill_r1_trimmed-corrected/contigs.fa outputs/glimmer/second_run.icm \
    #    outputs/glimmer/second_run.out \
    #    2>outputs/glimmer/second_run.err 1>&2
    return os.path.join(output_destination, '%s.glimmer.predict' % phage_name)


def run_trnascan_se(fasta_path, phage_name, output_destination):
    '''
    Run tRNAscan-SE and return output file
    '''
    print('*** RUNNING tRNAscan ***')

    out_file_path = os.path.join(output_destination, '%s.trnascan_se.results.txt' % phage_name)
    subprocess.run(['tRNAscan-SE', '-B', '-D', '-o', out_file_path, fasta_path], check=True)
    return out_file_path


# I think I want to add aragorn here in the hopes that it will be a
# bit better at picking up mtRNAs, which I imagine are useful for
# phages once the bacteria 'realize' there are shenanigans afoot and
# shut down translation.  In my brief reading of trnascan-se, I did
# not get the sense that it looks for them as effectively.
