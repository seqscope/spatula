# Use Ubuntu 22.04 as the base
FROM ubuntu:22.04

# Avoid interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install all required build tools and libraries
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    wget \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    autoconf \
    automake \
    libtool \
    libssl-dev \
    libcurl4-openssl-dev \
    libdeflate-dev \
    imagemagick \
 && rm -rf /var/lib/apt/lists/*

# Set the working directory
WORKDIR /app

# Copy the entire project into the container
COPY . .

# -------------------------------
# Build submodule: htslib
# -------------------------------
RUN cd submodules/htslib && \
    autoreconf -i && \
    ./configure && \
    make -j$(nproc)

# -------------------------------
# Build submodule: qgenlib
# -------------------------------
RUN cd submodules/qgenlib && \
    mkdir -p build && cd build && \
    cmake .. && \
    make -j$(nproc)

# -------------------------------
# Build main project: spatula
# -------------------------------
RUN mkdir -p build && cd build && \
    cmake .. && \
    make -j$(nproc)

# # Set default command to interactive shell
# CMD ["bash"]

SHELL ["/bin/bash", "-c"]

# Command to run when starting the container
COPY ./entrypoint.sh /
RUN chmod 755 /entrypoint.sh
ENTRYPOINT ["/entrypoint.sh"]
# Optionally set a default command (fallback if no args passed)
CMD []