FROM mambaorg/micromamba
LABEL maintainer="Anthony Cesnik <cesnik@wisc.edu>"

# Set up micromamba and base env
COPY --chown=root workflow/envs/spritzbase.yaml /tmp/spritzbase.yaml
RUN micromamba install -y -n base -f /tmp/spritzbase.yaml \
        && micromamba clean --all --yes

# install prebuilt SpritzCMD
COPY ./SpritzCMD/bin/x64/Release/net6.0/ /app/spritz/
WORKDIR /app/spritz/

