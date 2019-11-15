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

# Activating environment
# When running spritz, "conda activate spritz" should preceed a command
# I've had trouble running `conda activate` after `conda init` in the docker build
