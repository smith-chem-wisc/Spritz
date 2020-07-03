FROM conda/miniconda3
LABEL maintainer="Anthony Cesnik <cesnik@wisc.edu>"

# install dotnet
WORKDIR /usr/bin/local
RUN apt-get update -y && \
        apt-get install -y wget curl gpg \
        && wget -qO- https://packages.microsoft.com/keys/microsoft.asc | gpg --dearmor > microsoft.asc.gpg \
        && mv microsoft.asc.gpg /etc/apt/trusted.gpg.d/ \
        && wget -q https://packages.microsoft.com/config/debian/9/prod.list \
        && mv prod.list /etc/apt/sources.list.d/microsoft-prod.list \
        && chown root:root /etc/apt/trusted.gpg.d/microsoft.asc.gpg \
        && chown root:root /etc/apt/sources.list.d/microsoft-prod.list \
        && apt-get install -y apt-transport-https \
        && apt-get update -y \
        && apt-get install -y dotnet-sdk-2.2

# install Spritz
WORKDIR /app
COPY . ./
RUN conda init \
        && conda update conda \
        && conda env create --name spritz --file environment.yaml

# Add script to execute snakemake upon activate
# Note that running this script outside of activation has not worked
ADD start.sh /usr/local/envs/spritz/etc/conda/activate.d/start.sh
RUN chmod 777 /usr/local/envs/spritz/etc/conda/activate.d/start.sh

# activate environment
RUN echo "source activate spritz" > ~/.bashrc
ENV PATH /opt/conda/envs/spritz/bin:$PATH
