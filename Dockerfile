FROM bioconda/bioconda-utils-build-env
COPY main /usr/local/bin/phglab/ALFATClust
WORKDIR /home/work
RUN conda update -n base conda && \
    conda install -c conda-forge biopython python-igraph leidenalg numpy pandas scipy && \
    conda install -c bioconda mash mmseqs2 && \
    ln -s "/usr/local/bin/phglab/ALFATClust/alfatclust.py" "/usr/local/bin/alfatclust"
CMD [ "/bin/bash" ]
