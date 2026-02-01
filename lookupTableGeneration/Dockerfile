FROM ubuntu:24.04

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    tcl \
    tcl-dev \
    tk \
    tk-dev \
    libsuitesparse-dev \
    libplplot-dev \
    plplot-tcl \
    libopenblas-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy FLOOXS source (user must provide this folder)
COPY FLOOXS_2026/ /opt/flooxs/

# Set environment
ENV FLXSHOME=/opt/flooxs
ENV PL_LIBRARY=/usr/share/plplot5.15.0/tcl/

# Build FLOOXS
WORKDIR /opt/flooxs/src
RUN mkdir -p release && cd release && cmake .. && make -j$(nproc)

# Create working directory for simulations
WORKDIR /work

# The built executable
ENTRYPOINT ["/opt/flooxs/src/release/flooxs"]
