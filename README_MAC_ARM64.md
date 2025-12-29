# CGCT Installation Guide - macOS ARM64 (Apple Silicon M1/M2/M3)

This guide covers installing CGCT on macOS with Apple Silicon (ARM64) processors.

## Prerequisites

- macOS 11 (Big Sur) or later
- Apple Silicon Mac (M1, M2, M3, etc.)
- Administrator access for installing software

## Step 1: Install BLAST+

BLAST+ is required for sequence alignment and homology detection.

### Download BLAST+

1. Visit: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
2. For macOS ARM64, you have two options:
   - **`ncbi-blast-2.17.0+-aarch64-macosx.tar.gz`** (Recommended - latest)
   - **`ncbi-blast-2.17.0+-aarch64.dmg`** (Installer version)

### Install BLAST+ (tar.gz version - Recommended)

1. Download the `.tar.gz` file
2. Open Terminal
3. Navigate to your Downloads folder:
   ```bash
   cd ~/Downloads
   ```
4. Extract the archive:
   ```bash
   tar -xzf ncbi-blast-2.17.0+-aarch64-macosx.tar.gz
   ```
5. Move it to a standard location:
   ```bash
   sudo mv ncbi-blast-2.17.0+ /usr/local/ncbi-blast
   ```
6. Add it to your PATH by editing `~/.zprofile`:
   ```bash
   echo 'export PATH="/usr/local/ncbi-blast/bin:$PATH"' >> ~/.zprofile
   source ~/.zprofile
   ```

### Install BLAST+ (dmg version - Alternative)

1. Download the `.dmg` file
2. Double-click to mount it
3. Run the installer inside
4. Follow the installation wizard

### Verify Installation

Open Terminal and run:
```bash
blastn -version
```

You should see output showing the BLAST version.

## Step 2: Install Docker Desktop

Docker is required for SibeliaZ on macOS (SibeliaZ only natively runs on Linux).

### Download Docker Desktop

1. Visit: https://www.docker.com/products/docker-desktop
2. Click "Download for Mac"
3. Choose the **ARM64** version (Apple Silicon)
   - File will be named something like: `Docker.dmg` with "Apple Silicon" indicated

### Install Docker Desktop

1. Open the downloaded `.dmg` file
2. Drag the Docker icon to Applications
3. Open Applications and launch Docker.app
4. Follow the setup wizard
5. Enter your Mac password when prompted (Docker needs admin access)
6. Wait for Docker to fully start (you'll see the Docker menu icon in the top-right)

### Verify Installation

Open Terminal and run:
```bash
docker --version
```

You should see the Docker version number.

## Step 3: Set Up SibeliaZ Docker Image

SibeliaZ runs inside a Docker container on macOS.

### Option A: Pull Pre-built Docker Image (Easiest - Recommended)

1. Open Terminal
2. Pull the Docker image:
   ```bash
   docker pull aaronmauve/sibeliaz-conda:latest
   ```
3. This downloads the image (~2-3 GB) - takes a few minutes

### Option B: Build Docker Image Yourself

1. Make sure Docker Desktop is running
2. Get the `Dockerfile` from the CGCT distribution
3. Open Terminal in the directory containing the `Dockerfile`
4. Build the image:
   ```bash
   docker build -t sibeliaz-conda:latest .
   ```
   This can take several minutes depending on your internet connection

## Step 4: Configure CGCT

Edit the `config.txt` file (in the same directory as the CGCT program):

```
# These should be set by default, but verify:
USE_DOCKER_SIBELIAZ=false
DOCKER_SIBELIAZ_IMAGE=aaronmauve/sibeliaz-conda:latest
```

Save the file.

## Running CGCT

1. Make sure Docker Desktop is running (look for the Docker menu icon in the top-right)
2. Launch the CGCT program
3. The program will use:
   - BLAST+ from `/usr/local/ncbi-blast/bin` (or wherever you installed it)
   - SibeliaZ via Docker Desktop automatically
4. Load your sequences and run the analysis

## Troubleshooting

### "blastn not found"
- Check installation: `blastn -version` in Terminal
- If not found, verify the PATH:
  ```bash
  echo $PATH
  ```
- Should include `/usr/local/ncbi-blast/bin`
- If not, check that `~/.zprofile` contains the export line and you ran `source ~/.zprofile`

### "Docker not found" or "Cannot connect to Docker daemon"
- Make sure Docker Desktop is running
  - Check the top-right menu bar for the Docker icon
  - If not there, launch Docker.app from Applications
- Verify: `docker --version` in Terminal
- Docker may take a minute to start - wait and try again

### "docker pull" fails or times out
- Check internet connection
- Try again - Docker Hub sometimes has rate limits
- Verify Docker is running and has internet access

### "SibeliaZ timeout" or "exec format error"
- The ARM64 Docker image should handle this, but if you get errors:
  1. Verify the image exists: `docker images | grep sibeliaz`
  2. Try pulling the image again: `docker pull aaronmauve/sibeliaz-conda:latest`
  3. Check Docker Desktop has sufficient memory (Preferences → Resources)

### Program runs but no output
- Verify both BLAST+ and Docker are installed and working
- Check `config.txt` settings
- Try running from Terminal to see detailed error messages

## File Versions

The recommended versions above (BLAST+ 2.17.0) are current as of this guide. If newer versions are available at the NCBI FTP site, download the ARM64 version (`aarch64-macosx.tar.gz` or `aarch64.dmg`) instead.

## Performance Notes

- ARM64 Macs run the Docker image natively with good performance
- First run may take longer as Docker downloads the SibeliaZ image
- Subsequent runs will be faster (image is cached locally)

## Getting Help

If you encounter issues:
1. Verify installations: `blastn -version` and `docker --version` in Terminal
2. Ensure Docker Desktop is running (check menu bar)
3. Check that the Docker image exists: `docker images | grep sibeliaz`
4. Verify `config.txt` is in the right location and has correct settings
5. Check the program's debug output for specific error messages
