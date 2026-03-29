FROM condaforge/mambaforge:latest

# Install meme
RUN mamba install -y meme || echo "meme install failed, trying without"

# Create working directory
WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

EXPOSE 8001

CMD ["python", "server.py"]
