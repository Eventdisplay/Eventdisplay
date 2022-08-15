# Docker files and images for Eventdisplay

Docker images are made available for the following use cases.

(for the time being, there is some overlap with the [Eventdisplay container](https://github.com/Eventdisplay/Eventdisplay_Docker) repository.

## CTA prod5 analysis: analysis of prod5 CTA simulations

- [Dockerfile](dockerfiles/Dockerfile-cta-prod5) 
- docker image available from [container page](https://github.com/Eventdisplay/Eventdisplay/pkgs/container/eventdisplay) with tag `cta-prod5`

### Running

Sim_telarray files and output evndisp root files are read and writting from the ./data directory.

To analyse a prod5 sim_telarray file (replace VERSION by requested version):

```
$  docker run --rm -it -v "$(pwd)/data:/data" ghcr.io/eventdisplay/eventdisplay:VERSION-cta-prod5 \
     /workdir/run.sh \
     /data/gamma_20deg_0deg_run9___cta-prod5-paranal_desert-2147m-Paranal-dark_cone10.simtel.zst \
     CTA.prod5S.BL-4LSTs25MSTs70SSTs-MSTF.lis
```

To run the container in bash and analyse a prod5 sim_telarray file:

```
$ docker run --rm -it -v "$(pwd)/data:/data" ghcr.io/eventdisplay/eventdisplay:VERSION-cta-prod5 bash
$ run.sh \
    /data/gamma_20deg_0deg_run9___cta-prod5-paranal_desert-2147m-Paranal-dark_cone10.simtel.zst \
    CTA.prod5S.BL-4LSTs25MSTs70SSTs-MSTF.lis
```

### Building

Building expects that a tar ball of hessioxx (hessioxxx.tar.gz) is available in the building directory (GitHub action pulls it from a cloud directory).

```
$ docker build -f dockerfiles/Dockerfile-cta-prod5 -t eventdisplay-cta-dl1-prod5 .
```
(for rebuilt with `--no-cache`)

Download hessioxx from https://www.mpi-hd.mpg.de/hfm/CTA/MC/Software/Testing/hessioxxx.tar.gz (passwd applies)


## CTA slib: analysis library used for CTA

- [Dockerfile](dockerfiles/Dockerfile-cta-slib) 
- docker image available from [container page](https://github.com/Eventdisplay/Eventdisplay/pkgs/container/eventdisplay) with tag `cta-slib`

### Running

To run the container in bash (replace VERSION by requested version):

```
% docker run --rm -it -v "$(pwd)/:/workdir" ghcr.io/eventdisplay/eventdisplay:VERSION-cta-slib bash
```

### Building

```
% docker build -f dockerfiles/Dockerfile-cta-slib -t eventdisplay-cta-slib .
```

