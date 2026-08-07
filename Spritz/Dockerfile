# Pinned deliberately. This was previously untagged, so the image - and with it the micromamba and
# solver versions - drifted on every rebuild. That is how the build came to fail with no change to
# this repository. Bump this tag intentionally rather than inheriting a new base silently.
FROM mambaorg/micromamba:2.9.0-ubuntu24.04
LABEL maintainer="Anthony Cesnik <cesnik@wisc.edu>"

# Set up micromamba and base env
COPY --chown=root workflow/envs/spritzbase.yaml /tmp/spritzbase.yaml
RUN micromamba install -y -n base -f /tmp/spritzbase.yaml \
        && micromamba clean --all --yes

# install prebuilt SpritzCMD
COPY ./SpritzCMD/bin/x64/Release/net10.0/ /app/spritz/
WORKDIR /app/spritz/

