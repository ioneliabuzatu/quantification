#! /usr/bin/env python3
import os
import sys
import subprocess
from Bio import SeqIO
from multiprocessing import Pool, cpu_count
from glob import glob
import shutil
import plot as plot
import matplotlib.pyplot as plt



def detect(seq_dir, rootdir, database, barcode, pathogen):


    max_threads = cpu_count()
    free = True
    abs_rootdir = os.path.abspath(rootdir)
    seq_rootdir = os.path.abspath(seq_dir)
    batch = 0
    files_old = []
    while free:
        batch += 1
        fastq = [y for x in os.walk(seq_rootdir) for y in glob(os.path.join(x[0], '*.fastq')) if "analysed" not in y]
        fq = [y for x in os.walk(seq_rootdir) for y in glob(os.path.join(x[0], '*.fq')) if "analysed" not in y]
        files = fastq + fq
        file_to_do = []
        for file in files:
            if file not in files_old:
                file_to_do.append(file)
        files_old = files_old + file_to_do

        output_barcode = os.path.join(abs_rootdir, "analysed_barcode")
        analysed_fasta = os.path.join(abs_rootdir, "analysed_fasta")
        analysed_fastq = os.path.join(abs_rootdir, "analysed_fastq")
        analysed_blastn = os.path.join(abs_rootdir, "analysed_blastn")
        analysing_fasta = os.path.join(abs_rootdir, "analysing_fasta")
        if not os.path.exists(output_barcode):
            os.mkdir(output_barcode)
        if not os.path.exists(analysed_fasta):
            os.mkdir(analysed_fasta)
        if not os.path.exists(analysed_fastq):
            os.mkdir(analysed_fastq)
        if not os.path.exists(analysed_blastn):
            os.mkdir(analysed_blastn)
        if not os.path.exists(analysing_fasta):
            os.mkdir(analysing_fasta)

        if len(file_to_do) == 0:
            sys.exit("NO NEW FASTQ FILES")
        fastq_combined = os.path.join(analysed_fastq, 'combined.' + str(batch) + '.fq')
        with open(fastq_combined, 'wb') as wfd:
            for fileIn in file_to_do:
                if fileIn.endswith("fastq") or fileIn.endswith("fq"):
                    if not "combined" in fileIn:
                        with open(fileIn, 'rb') as fd:
                            shutil.copyfileobj(fd, wfd, 1024 * 1024 * 10)
        sys.stdout.write("ANALYSING BATCH N. " + str(batch) + "\n")

        #for file in file_to_do:
        #    shutil.move(file, os.path.join(analysed_fastq, file.split("/")[-1]))

        if barcode:
            with open(os.path.join(abs_rootdir, "porechop.log"), "w") as fhlog, open(
                    os.path.join(abs_rootdir, "porechop.err"), "w") as fherr:
                cmd_porechop = "porechop -i %s -b %s" % (fastq_combined, output_barcode)
                p = subprocess.Popen(cmd_porechop, stdout=fhlog, stderr=fherr, universal_newlines=True, shell=True)
                p.communicate()
        reads_number = []
        for root, dirs, files in os.walk(output_barcode):
            for file in files:
                if file.startswith("BC") and file.endswith("fastq"):
                    count = 0
                    file_split = os.path.join(root, file)
                    fasta_name = ".".join([file.split(".")[0], str(batch), "fasta"])
                    fasta_seq = []
                    for record in SeqIO.parse(file_split, "fastq"):
                        count += 1
                        record.id = file.split(".")[0] + "_" + str(count)
                        fasta_seq.append(record)
                    reads_number.append([file.split(".")[0], count])
                    SeqIO.write(fasta_seq, os.path.join(analysing_fasta, fasta_name), "fasta")

        fasta = [y for x in os.walk(analysing_fasta) for y in glob(os.path.join(x[0], '*.fasta'))]
        fasta_ready = []
        for fasta_file in fasta:
            fasta_ready.append([fasta_file, database, analysed_blastn])

        if len(fasta) < max_threads:
            max_threads = len(fasta)
        with Pool(max_threads) as p:
            results = p.map(blast, fasta_ready)
        for file_fasta in fasta:
            shutil.move(file_fasta, os.path.join(analysed_fasta, file_fasta.split("/")[-1]))

        plot.barplot(results)#, reads_number)
        plt.show()

    print("DONE")


def blast(data):
    blast_out_name = data[0].split("/")[-1].split(".")[0] + ".blastn"
    blast_out_path = os.path.join(data[2], blast_out_name)
    cmd_bast = "blastn -db %s -query %s -max_target_seqs 1 -outfmt  \" 6 qseqid ssciname \" -out /dev/stdout " % (data[1], data[0])
    p = subprocess.Popen(cmd_bast, stdout=subprocess.PIPE, shell=True)
    out = "\n".join(list(set(p.communicate()[0].decode().split("\n"))))
    if os.path.exists(blast_out_path):
        append_write = 'a'
    else:
        append_write = 'w'
    with open(blast_out_path, append_write) as fh:
        fh.write(out + "\n")
    return blast_out_path


if __name__ == '__main__':
    detect(*sys.argv[1:])