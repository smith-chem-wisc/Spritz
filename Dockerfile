FROM ubuntu:latest
LABEL maintainer="Anthony Cesnik <cesnik@wisc.edu>"

ENV TZ=America/Aruba

# install dotnet core
WORKDIR /usr/bin/local
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone && \
	apt-get update -y && \
	apt-get install -y wget curl software-properties-common sudo && \
	useradd -m docker && echo "docker:docker" | chpasswd && adduser docker sudo && \
	wget -q https://packages.microsoft.com/config/ubuntu/18.04/packages-microsoft-prod.deb && \
	dpkg -i packages-microsoft-prod.deb && \
	add-apt-repository universe && \
	apt-get install -y apt-transport-https && \
	apt-get update -y && \
	apt-get install -y dotnet-sdk-2.2

# install Spritz
WORKDIR /app
COPY . ./
RUN dotnet restore && \
	dotnet build CMD/CMD.csproj && \
	dotnet /app/CMD/bin/Debug/netcoreapp2.1/CMD.dll -c setup
WORKDIR /app/CMD/bin/Debug/netcoreapp2.1