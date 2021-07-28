FROM mambaorg/micromamba:0.15.2
LABEL maintainer="Anthony Cesnik <cesnik@wisc.edu>"

# Set up micromamba and base env
COPY --chown=micromamba:micromamba workflow/envs/spritzbase.yaml /tmp/spritzbase.yaml
RUN micromamba install -y -n base -f /tmp/spritzbase.yaml \
        && micromamba clean --all --yes

# install Spritz
WORKDIR /app/spritz/workflow
COPY . ../

# activate environment
RUN echo "source activate" > ~/.bashrc
