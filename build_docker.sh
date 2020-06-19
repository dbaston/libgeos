#!/usr/bin/env bash
set -euxo pipefail

docker build -t dbaston/geos-shapely:latest -f Dockerfile.pytest .
docker push  dbaston/geos-shapely:latest

docker build -t dbaston/geos-postgis:latest -f Dockerfile.pgtest .
docker push dbaston/geos-postgis:latest

docker build -t dbaston/geos-sf:latest -f Dockerfile.rtest .
docker push dbaston/geos-sf:latest
