FROM debian:latest

RUN apt update && apt upgrade -y
RUN apt install curl unzip -y

# Download Atria and install its dependencies
RUN curl -L -O https://github.com/cihga39871/Atria/releases/download/v4.0.3/atria-4.0.3-linux-ubuntu22.tar.gz
RUN tar -zxf atria-4.0.3-linux-ubuntu22.tar.gz
RUN rm atria-4.0.3-linux-ubuntu22.tar.gz
RUN apt install pigz pbzip2 -y
ENV PATH="/atria-4.0.3/bin:$PATH"

# Download FASTQC and install its dependencies (JRE and Perl)
RUN curl -L -O https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip
RUN unzip fastqc_v0.12.1.zip
RUN rm fastqc_v0.12.1.zip
RUN apt install default-jre perl -y

# Download featureCounts
RUN curl -L -O https://sourceforge.net/projects/subread/files/subread-2.0.6/subread-2.0.6-Linux-x86_64.tar.gz
RUN tar -zxf /subread-2.0.6-Linux-x86_64.tar.gz
RUN rm subread-2.0.6-Linux-x86_64.tar.gz
ENV PATH="/subread-2.0.6-Linux-x86_64/bin:$PATH"

# Download STAR and compile it (also install required dependencies)
RUN curl -L -O https://github.com/alexdobin/STAR/archive/2.7.11b.tar.gz
RUN tar -xzf 2.7.11b.tar.gz
RUN rm 2.7.11b.tar.gz
RUN apt-get install build-essential vim-common zlib1g-dev xxd -y
RUN make -C /STAR-2.7.11b/source STAR
ENV PATH="/STAR-2.7.11b/bin/Linux_x86_64:$PATH"

# Download and compile samtools, also install dependencies
RUN curl -L -O https://github.com/samtools/samtools/releases/download/1.20/samtools-1.20.tar.bz2
RUN tar -xvjf samtools-1.20.tar.bz2
RUN rm samtools-1.20.tar.bz2
RUN apt-get install libncurses-dev liblzma-dev  libbz2-dev libcurl4-openssl-dev -y
RUN make -C /samtools-1.20
RUN make -C /samtools-1.20 install

# Install software-properties-common, which also includes Python
RUN apt-get install software-properties-common -y
RUN apt-get update

# Install bedtools
RUN curl -L -O https://github.com/arq5x/bedtools2/archive/refs/tags/v2.31.1.zip
RUN unzip v2.31.1.zip
RUN rm v2.31.1.zip
RUN make -C /bedtools2-2.31.1
ENV PATH="/bedtools2-2.31.1/bin:$PATH"



