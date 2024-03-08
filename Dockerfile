#
# Docker image to run sequence to taxonomy analysis
#
# Includes reference data (blastdb)
#
# NEEDS
#   python3
#   Anaconda ? 
#     pandas
#     blast

# Ubunutu base
FROM ubuntu:22.04

# make sure installs dont hang on user input
ENV DEBIAN_FRONTEND noninteractive


# install r-base and pre-requisitis
# installs R 3.6.3 (on ubuntu:20.04)
RUN set -e \
      && apt-get -y update \
      && apt-get -y dist-upgrade \
      && apt-get -y install --no-install-recommends --no-install-suggests \
        apt-transport-https apt-utils ca-certificates locales pandoc pkg-config \
        ncbi-blast+ \
      && apt-get -y autoremove \
      && apt-get clean \
      && rm -rf /var/lib/apt/lists/*

# install several Python packages with Conda p
RUN mkdir ./conda
RUN ./create_conda_env.sh

# UTF-8 mode
RUN set -e \
      && locale-gen en_US.UTF-8 \
      && update-locale

#
# copy in our application
#
# do this as a git clone, instead!?!?
COPY classify_sequence .
COPY version_git.txt .

#
# copy in reference data
#
RUN mkdir -p ./blast
COPY blast/ ./blast/
RUN find ./blast/
COPY processed_accessions_e.tsv ./
COPY processed_accessions_a.tsv ./
COPY fixed_vmr_a.tsv ./
COPY fixed_vmr_e.tsv ./

# what does ENTRYPOINT do exactly?
CMD [ "./classify_sequence" ]
