FROM buisciii/centos7_base_image:latest

COPY ./scif_app_recipes/* /opt/

RUN echo "Install basic development tools" && \
    yum -y groupinstall "Development Tools" && \
    yum -y update && yum -y install wget curl && \
    echo "Install python2.7 setuptools and pip" && \
    yum -y install python-setuptools && \
    easy_install pip && \
    echo "Installing SCI-F" && \
    pip install scif ipython

RUN echo "Installing FastQC app" && \
    scif install /opt/fastqc_v0.11.7_centos7.scif && \
    echo "Installing trimmomatic app" && \
    scif install /opt/trimmomatic_v0.38_centos7.scif && \
    echo "Installing samtools app" && \
    scif install /opt/samtools_v1.9_centos7.scif && \
    echo "Installing bedtools app" && \
    scif install /opt/bedtools_v2.27_centos7.scif && \
    echo "Installing varscan app" && \
    scif install /opt/varscan_v2.3.9_centos7.scif && \
    echo "Installing multiqc app" && \
    scif install /opt/multiqc_v1.4_centos7.scif && \
    echo "Installing bwa app" && \
    scif install /opt/bwa_v0.7.17_centos7.scif && \
    echo "Installing kggseq app" && \
    scif install /opt/kggseq_v1.1_centos7.scif

    ## R packages

    # Install core R dependencies
    RUN echo "r <- getOption('repos'); r['CRAN'] <- 'https://ftp.acc.umu.se/mirror/CRAN/'; options(repos = r);" > ~/.Rprofile

# Include ENV variables
ENV LC_ALL=en_US.UTF-8
ENV PATH=$PATH:/scif/apps/aragorn/bin
ENV PATH=$PATH:/scif/apps/barrnap/bin
ENV PATH=$PATH:/scif/apps/bedtools/bin
ENV PATH=$PATH:/scif/apps/bowtie2/bin
ENV PATH=$PATH:/scif/apps/bwa/bin
ENV PATH=$PATH:/scif/apps/fastqc/bin
ENV PATH=$PATH:/scif/apps/gcc/bin
ENV PATH=$PATH:/scif/apps/hmmer3/bin
ENV PATH=$PATH:/scif/apps/htslib/bin
ENV PATH=$PATH:/scif/apps/kggseq/bin
ENV PATH=$PATH:/scif/apps/minced/bin
ENV PATH=$PATH:/scif/apps/multiqc/bin
ENV PATH=$PATH:/scif/apps/ncbiblast/bin
ENV PATH=$PATH:/scif/apps/picard/bin
ENV PATH=$PATH:/scif/apps/pilon/bin
ENV PATH=$PATH:/scif/apps/prodigal/bin
ENV PATH=$PATH:/scif/apps/prokka/bin
ENV PATH=$PATH:/scif/apps/python3/bin
ENV PATH=$PATH:/scif/apps/quast/bin
ENV PATH=$PATH:/scif/apps/samtools/bin
ENV PATH=$PATH:/scif/apps/spades/bin
ENV PATH=$PATH:/scif/apps/sratoolkit/bin
ENV PATH=$PATH:/scif/apps/tbl2asn/bin
ENV PATH=$PATH:/scif/apps/trimmomatic/bin
ENV PATH=$PATH:/scif/apps/unicycler/bin
ENV PATH=$PATH:/scif/apps/varscan/bin
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scif/apps/aragorn/lib/lib
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scif/apps/barrnap/lib/lib
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scif/apps/bedtools/lib/lib
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scif/apps/bowtie2/lib/lib
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scif/apps/bwa/lib/lib
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scif/apps/fastqc/lib/lib
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scif/apps/gcc/lib/lib
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scif/apps/hmmer3/lib/lib
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scif/apps/htslib/lib/lib
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scif/apps/kggseq/lib/lib
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scif/apps/minced/lib/lib
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scif/apps/multiqc/lib/lib
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scif/apps/ncbiblast/lib/lib
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scif/apps/picard/lib/lib
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scif/apps/pilon/lib/lib
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scif/apps/prodigal/lib/lib
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scif/apps/prokka/lib/lib
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scif/apps/python3/lib/lib
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scif/apps/quast/lib/lib
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scif/apps/samtools/lib/lib
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scif/apps/spades/lib/lib
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scif/apps/sratoolkit/lib/lib
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scif/apps/tbl2asn/lib/lib
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scif/apps/trimmomatic/lib/lib
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scif/apps/unicycler/lib/lib
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scif/apps/varscan/lib/lib
#ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH/usr/local/lib

#ENTRYPOINT ["/opt/docker-entrypoint.sh"]
#CMD ["scif"]
RUN echo "export LC_ALL=en_US.UTF-8" >> /etc/bashrc
RUN find /scif/apps -maxdepth 2 -name "bin" | while read in; do echo "export PATH=\$PATH:$in" >> /etc/bashrc;done
RUN if [ -z "${LD_LIBRARY_PATH-}" ]; then echo "export LD_LIBRARY_PATH=/usr/local/lib" >> /etc/bashrc;fi
RUN find /scif/apps -maxdepth 2 -name "lib" | while read in; do echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:$in" >> /etc/bashrc;done
