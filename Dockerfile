FROM mambaorg/micromamba:2.3.0

USER root

WORKDIR /opt/PKS_DOCK

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml requirements-pip.lock /tmp/

RUN micromamba install -y -n base -f /tmp/environment.yml \
    && micromamba clean --all --yes \
    && rm -f /tmp/environment.yml /tmp/requirements-pip.lock

ARG MAMBA_DOCKERFILE_ACTIVATE=1

COPY --chown=$MAMBA_USER:$MAMBA_USER . /opt/PKS_DOCK

RUN download-antismash-databases

RUN python -m pip install --no-deps .

RUN mkdir -p /opt/PKS_DOCK/results \
    /opt/PKS_DOCK/logs \
    /opt/PKS_DOCK/cache

WORKDIR /opt/PKS_DOCK

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "pks-dock"]
