# Docker build context should be one directory up from the "docker" directory. I.e.
#
#     docker build -t gmxapi/docs -f docs.dockerfile ..
#
# Note that this image depends on the image gmxapi/ci-mpich, which is built from
# ci.dockerfile. Build errors will occur if the current repository is too
# different from the version used to build gmxapi/ci-mpich. Either (re)build a
# local copy of gmxapi/ci-mpich or specify a particular tagged version
# gmxapi/ci-mpich:$REF by passing the build argument REF set to one of the tags
# at https://cloud.docker.com/u/gmxapi/repository/docker/gmxapi/ci-mpich
# I.e.
#
#     docker build -t gmxapi/docs -f docs.dockerfile --build-arg REF=sometag ..
#
# Launch and bind the host port 8080 to the web server in the container.
#
#     docker run --rm -p 8080:80 gmxapi/docs
#
# Connect by browsing to http://localhost:8080/
#

ARG REF=latest
FROM gmxapi/ci-mpich:$REF as docsbuild

USER root

RUN apt-get update && \
    apt-get -yq --no-install-suggests --no-install-recommends install \
        plantuml && \
    rm -rf /var/lib/apt/lists/*

USER testing

RUN . $VENV/bin/activate && \
    pip install -r /home/testing/gmxapi/requirements-docs.txt --no-cache-dir

COPY documentation /home/testing/gmxapi/documentation
RUN cd /home/testing/gmxapi && \
    . $VENV/bin/activate && \
    sphinx-build -b html documentation html


FROM httpd

COPY --from=docsbuild /home/testing/gmxapi/html/ /usr/local/apache2/htdocs/
