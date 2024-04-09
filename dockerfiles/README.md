# Docker files and images for Eventdisplay

Docker images are made available for the following use cases.

## CTA prod5/prod6 analysis analysis: analysis of prod5 CTA simulations

(for prod5, just replace prod6 in the following in all steps)

- [Dockerfile](dockerfiles/Dockerfile-cta-prod6)
- docker image available from [container page](https://github.com/Eventdisplay/Eventdisplay/pkgs/container/eventdisplay) with tag `cta-prod6`

### Running prod5/prod6

Sim_telarray files and output evndisp root files are read and writing from the ./data directory.

To analyze a prod6 sim_telarray file (replace VERSION by requested version):

```bash
$  docker run --rm -it -v "$(pwd)/data:/data" ghcr.io/eventdisplay/eventdisplay:VERSION-cta-prod6 \
     /home/run.sh \
     /data/gamma_20deg_0deg_run9___cta-prod6-paranal_desert-2147m-Paranal-dark_cone10.simtel.zst
```

To run the container in bash and analyze a prod6 sim_telarray file:

```bash
$ docker run --rm -it -v "$(pwd)/data:/data" ghcr.io/eventdisplay/eventdisplay:VERSION-cta-prod6 bash
$ run.sh \
    /data/gamma_20deg_0deg_run9___cta-prod6-paranal_desert-2147m-Paranal-dark_cone10.simtel.zst
```

### Building prod5/prod6

Building expects that a tar ball of hessioxx (hessioxxx.tar.gz) is available in the building directory (GitHub action pulls it from a cloud directory).

```bash
docker build -f dockerfiles/Dockerfile-cta-prod6 -t eventdisplay-cta-dl1-prod6 .
```

(for rebuilt with `--no-cache`)

Download hessioxx from https://www.mpi-hd.mpg.de/hfm/CTA/MC/Software/Testing/hessioxxx.tar.gz (passwd applies)

## CTA slib: analysis library used for CTA

- [Dockerfile](dockerfiles/Dockerfile-cta-slib)
- docker image available from [container page](https://github.com/Eventdisplay/Eventdisplay/pkgs/container/eventdisplay) with tag `cta-slib`

### Running cta-slib

To run the container in bash (replace VERSION by requested version):

```bash
docker run --rm -it -v "$(pwd)/:/workdir" ghcr.io/eventdisplay/eventdisplay:VERSION-cta-slib bash
```

### Building cta-slib

```bash
% docker build -f dockerfiles/Dockerfile-cta-slib -t eventdisplay-cta-slib .
```
