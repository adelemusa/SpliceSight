FROM continuumio/miniconda3:4.12.0

# Add bioconda channels
RUN conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    conda config --set channel_priority strict

# Install meme with specific python version
RUN conda create -y -n meme_env python=3.9 && \
    conda clean -afy

# Activate environment and install meme
SHELL ["conda", "run", "-n", "meme_env", "/bin/bash", "-c"]
RUN conda install -y meme && conda clean -afy

# Verify meme installation
RUN conda run -n meme_env which fimo

# Create working directory
WORKDIR /app

COPY requirements.txt .
RUN conda run -n meme_env pip install --no-cache-dir -r requirements.txt

COPY . .

EXPOSE 8001

# Set up environment
ENV PATH="/opt/conda/envs/meme_env/bin:${PATH}"
ENV PYTHONPATH="/app"

CMD ["python", "server.py"]
