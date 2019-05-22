Bootstrap: docker
From: buisciii/centos7_base_image:latest

%files
    ./scif_app_recipes/ /opt/
%post
    echo "Install basic development tools"
    yum -y groupinstall "Development Tools"
    yum -y update && yum -y install wget curl openssl-devel geos-devel udunits2-devel libxml2-devel cairo-devel libgit2-devel


    echo "Install python2.7 setuptools and pip"
    yum -y install python-setuptools
    easy_install pip

    echo "Installing SCI-F"
    pip install scif

    echo "Installing FastQC app" && \
    scif install /opt/scif_app_recipes/fastqc_v0.11.7_centos7.scif && \
    echo "Installing trimmomatic app" && \
    scif install /opt/scif_app_recipes/trimmomatic_v0.38_centos7.scif && \
    echo "Installing bwa app" && \
    scif install /opt/scif_app_recipes/bwa_v0.7.17_centos7.scif && \
    echo "Installing picard app" && \
    scif install /opt/scif_app_recipes/picard_v1.140_centos7.scif && \
    echo "Installing htslib app" && \
    scif install /opt/scif_app_recipes/htslib_v1.9_centos7.scif && \
    echo "Installing bedtools app" && \
    scif install /opt/scif_app_recipes/bedtools_v2.27_centos7.scif && \
    echo "Installing samtools app" && \
    scif install /opt/scif_app_recipes/samtools_v1.9_centos7.scif && \
    echo "Installing varscan app" && \
    scif install /opt/scif_app_recipes/varscan_v2.3.9_centos7.scif && \
    echo "Installing bcftools app" && \
    scif install /opt/scif_app_recipes/bcftools_v1.9_centos7.scif && \
    echo "Installing kggseq app" && \
    scif install /opt/scif_app_recipes/kggseq_v1.1_centos7.scif && \
    echo "Installing R app" && \
    scif install /opt/scif_app_recipes/R_v3.5.1_centos7.scif && \
    echo "Installing multiqc app" && \
    scif install /opt/scif_app_recipes/multiqc_v1.8dev_centos7.scif && \
    echo "Installing bamutil app" && \
    scif install /opt/scif_app_recipes/bamutil_v1.0.13_centos7.scif
	


    # Executables must be exported for nextflow, if you use their singularity native integration.
    # It would be cool to use $SCIF_APPBIN_bwa variable, but it must be set after PATH variable, because I tried to use it here and in %environment without success.
    find /scif/apps -maxdepth 2 -name "bin" | while read in; do echo "export PATH=\${PATH}:$in" >> $SINGULARITY_ENVIRONMENT;done

    find /scif/apps -maxdepth 2 -name "lib" | while read in; do echo "export LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:$in" >> $SINGULARITY_ENVIRONMENT ;done

    if [[ ":$PATH:" == *":/scif/apps/snppipeline:"* ]];then

        export CLASSPATH=/scif/apps/varscan/varscan-2.3.9/varscan-2.3.9.jar:$CLASSPATH >> $SINGULARITY_ENVIRONMENT
        export CLASSPATH=/scif/apps/picard/picard.jar:$CLASSPATH >> $SINGULARITY_ENVIRONMENT

    fi

%runscript
    exec scif "$@"
