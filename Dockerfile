FROM bioconda/bioconda-utils-build-env
COPY main /usr/local/bin/phglab/ALFATClust
WORKDIR /home/work
RUN conda update -n base conda && \
    conda install -c conda-forge biopython=1.78 python-igraph=0.8.3 leidenalg=0.8.3 numpy=1.20.1 pandas=1.2.2 scipy=1.6.0 && \
    conda install -c bioconda mash=2.2.2 mmseqs2=12.113e3 && \
    ln -s "/usr/local/bin/phglab/ALFATClust/alfatclust.py" "/usr/local/bin/alfatclust"
CMD [ "/bin/bash" ]
