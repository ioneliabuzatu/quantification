FROM ubuntu:xenial

RUN apt-get clean all && apt-get update && apt-get install -y -q build-essential git wget perl \
    python3.5 python2.7 software-properties-common python3-pip sudo cmake samtools bedtools zlib1g-dev libc6 aptitude \
    libdbd-mysql-perl libdbi-perl libboost-all-dev libncurses5-dev bowtie default-jre parallel nano ncbi-blast+


RUN pip3 install biopython numpy  bcbio-gff==0.6.4 pandas==0.19.1 pybedtools==0.7.8 gffutils regex pysam matplotlib progressbar2 \
    psutil memory_profiler pathlib colorama simplesam tqdm Flask

WORKDIR /opt/

RUN git clone https://github.com/rrwick/Porechop.git && cd Porechop && python3 setup.py install

WORKDIR /opt/
