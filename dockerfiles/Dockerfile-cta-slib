From ubuntu:20.04

LABEL maintainer.name="Eventdisplay Team"
LABEL maintainer.email="gernot.maier@desy.de"

# force ubunto to use bash for /bin/sh
# (otherwise trouble with many evndisp scripts)
RUN yes no | dpkg-reconfigure dash

# Basic packages
RUN apt-get update && apt-get install -y \
  bash \
  bzip2 \
  gcc \
  g++ \
  git \
  libblas-dev \
  libbz2-dev \
  libfreetype6 \
  libxpm4 \
  libtiff-dev \
  libxtst6 \
  make \
  unzip \
  vim \
  wget \
  zstd
RUN apt-get clean

WORKDIR /home/
ENV EVNDISP /home/

# Install ROOT
ARG ROOTVERSION=6.24.06
RUN wget https://root.cern/download/root_v$ROOTVERSION.Linux-ubuntu20-x86_64-gcc9.3.tar.gz && \
    tar -xzf root_v$ROOTVERSION.Linux-ubuntu20-x86_64-gcc9.3.tar.gz && \
    rm -f root_v$ROOTVERSION.Linux-ubuntu20-x86_64-gcc9.3.tar.gz
ENV ROOTSYS ${EVNDISP}/root
ENV LD_LIBRARY_PATH "${ROOTSYS}/lib/:${LD_LIBRARY_PATH}"
ENV PATH ${ROOTSYS}/bin:${PATH}

# Install Eventdisplay
ARG EVNVERSION=5.6.7
ENV EVNDISPSYS ${EVNDISP}/Eventdisplay/
ENV SOFASYS ${EVNDISPSYS}/sofa/
RUN wget https://github.com/Eventdisplay/Eventdisplay/archive/refs/tags/v$EVNVERSION.zip && \
    unzip v$EVNVERSION.zip && rm -f v$EVNVERSION.zip && \
    mv Eventdisplay-$EVNVERSION Eventdisplay && \
    cd ${EVNDISPSYS} && ls -l && ./install_sofa.sh && \
    make -j 4 slib && strip -g lib/libVAnaSum.so &&  \
    rm -f obj/*.o obj/*.cpp

#RUN wget https://github.com/Eventdisplay/Eventdisplay/archive/refs/heads/root-cpp-versions.zip && \
#    unzip root-cpp-versions.zip && rm -f root-cpp-versions.zip && \
#    mv Eventdisplay-root-cpp-versions Eventdisplay && \
#    cd ${EVNDISPSYS} && ls -l && ./install_sofa.sh && \
#    make -j 2 slib && strip -g lib/libVAnaSum.so &&  \
#    rm -f obj/*.o obj/*.cpp

ENV LD_LIBRARY_PATH "${EVNDISPSYS}/obj:${LD_LIBRARY_PATH}"
ENV ROOT_INCLUDE_PATH ${EVNDISPSYS}/inc

# Install Eventdisplay Analysis files
ARG ANAVERSION=1.5.2
RUN wget https://github.com/Eventdisplay/Eventdisplay_AnalysisFiles_CTA/archive/refs/tags/v$ANAVERSION.zip && \
    unzip v$ANAVERSION.zip && rm -f v$ANAVERSION.zip && \
    mv Eventdisplay_AnalysisFiles_CTA-$ANAVERSION Eventdisplay_AnalysisFiles_CTA
ENV CTA_EVNDISP_AUX_DIR /home/Eventdisplay_AnalysisFiles_CTA/
