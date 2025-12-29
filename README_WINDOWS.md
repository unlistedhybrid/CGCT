# CGCT Installation Guide - Windows (x64)

This guide covers installing CGCT on Windows with all required dependencies.

## Prerequisites

- Windows 7 or later (64-bit)
- Administrator access for installing software

## Step 1: Install BLAST+

BLAST+ is required for sequence alignment and homology detection.

### Download BLAST+

1. Visit: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
2. Look for the file matching your system. For Windows 64-bit, download:
   - **`ncbi-blast-2.17.0+-win64.exe`** (or any newer version available)

### Install BLAST+

1. Download the `.exe` file
2. Double-click to run the installer
3. Follow the installation wizard (default options are fine)
4. Installation location: `C:\Program Files\NCBI\blast\bin` (default)

### Verify Installation

Open PowerShell and run:
```powershell
blastn -version
```

You should see output showing the BLAST version. If not, the `bin` folder may not be in your PATH - contact support.

## Step 2: Configure SibeliaZ - Choose Your Method

SibeliaZ (genome alignment tool) can run on Windows in two ways. Choose one:

If you have Windows 10 or newer the WSL option is simpler, it doesn't require downloading Docker software.
Running the WSL command for the first time during setup can take several minutes but it is simple and you are only required to do it one time.
If you have an older version of Windows or you are already familiar with Docker, you may prefer the Docker option.

If you don't know which option to pick, go with Option A with WSL and skip Option B.

---

### Option A: WSL + SibeliaZ (Recommended)

WSL (Windows Subsystem for Linux) allows Windows to run Linux tools natively. This is the recommended option for most users.

#### Option A1: Automatic Installation (Easiest)

This CGCT package includes an automatic installer script.

1. Make sure WSL is installed:
   - Open PowerShell as Administrator
   - Run:
     ```powershell
     wsl --install
     ```
   - Restart your computer if prompted

2. After extracting the CGCT folder, double-click:
   ```
   install-sibeliaz-wsl.bat
   ```

3. The installer will:
   - Launch WSL
   - Create a conda environment
   - Install SibeliaZ automatically

4. When finished, you should see a success message.

**Note:** This process may take several minutes the first time, but only needs to be done once.

---

#### Option A2: Manual Installation (If Automatic Fails)

WSL (Windows Subsystem for Linux) lets Windows run Linux tools natively. This is the **recommended option** for most users.

#### Install WSL2

1. Open PowerShell as Administrator
2. Run:
   ```powershell
   wsl --install
   ```
3. This installs WSL2 and Ubuntu by default
4. Restart your computer if prompted
5. After restart, Ubuntu will complete setup automatically

#### Install SibeliaZ in WSL

1. Open PowerShell and enter WSL:
   ```powershell
   wsl
   ```
   You're now in a Linux environment (notice the prompt changes)

2. Update package manager:
   ```bash
   sudo apt-get update
   ```

3. Install conda:
   ```bash
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   bash Miniconda3-latest-Linux-x86_64.sh
   ```
   Follow the installer prompts

4. Restart WSL (type `exit` then open PowerShell again and run `wsl`)

5. Install SibeliaZ:
   ```bash
   conda install -y -c bioconda sibeliaz twopaco
   ```

6. Test it:
   ```bash
   sibeliaz
   ```

7. If it is configured correctly it will show "You must provide the input file name"

---

### Option B: Docker + SibeliaZ (Alternative)

Docker containers provide consistent, isolated environments. Use this if you're already familiar with Docker or prefer Docker workflows.

#### Install Docker Desktop (if not already installed)

1. Visit: https://www.docker.com/products/docker-desktop
2. Click "Download for Windows"
3. Choose the installer for Windows (AMD64 is the standard for Windows)
4. Run the downloaded installer
5. Follow the installation wizard
6. When prompted, enable "WSL 2" (recommended)
7. Restart your computer when installation completes

#### Set Up SibeliaZ Docker Image

**Option B1: Pull Pre-built Image (Easiest)**

1. Open PowerShell
2. Pull the Docker image:
   ```powershell
   docker pull aaronmauve/sibeliaz-conda:latest
   ```
3. Test it:
   ```powershell
   docker run --rm aaronmauve/sibeliaz-conda:latest sibeliaz
   ```

4. If it is configured correctly it will show "You must provide the input file name"

**Option B2: Build Docker Image Yourself**

1. Get the `Dockerfile` from the CGCT distribution
2. Open PowerShell in the directory containing the `Dockerfile`
3. Build the image:
   ```powershell
   docker build -t sibeliaz-conda:latest .
   ```
   This takes ~5-10 minutes

## Step 3: Final Configuration

Edit `config.txt` with your chosen method:

**If you chose WSL (Option A):**
```
USE_DOCKER_SIBELIAZ=false
WSL_DISTRO=Ubuntu
WSL_USER=
```
The defaults work if you followed the installation instructions above. Only modify `WSL_DISTRO` or `WSL_USER` if you used a different Linux distribution or created a specific WSL user.

**If you chose Docker (Option B):**
```
USE_DOCKER_SIBELIAZ=true
DOCKER_SIBELIAZ_IMAGE=aaronmauve/sibeliaz-conda:latest
```
Change `DOCKER_SIBELIAZ_IMAGE` to `sibeliaz-conda:latest` if you built the Docker image yourself in Option B2 (or use whatever image name you chose).

Save the file.

## Running CGCT

Once all dependencies are installed and configured:

1. Launch the CGCT program
2. The program should automatically find:
   - BLAST+ executables
   - SibeliaZ (either via WSL or Docker, depending on your config)
   - progressiveMauve included with CGCT in the Tools folder
3. Load your sequences and run the analysis

## Troubleshooting

### "blastn not found"
- Check BLAST+ is installed: `blastn -version`
- If not found, add `C:\Program Files\NCBI\blast\bin` to your Windows PATH:
  1. Press `Win + X`, select "System"
  2. Click "Advanced system settings"
  3. Click "Environment Variables"
  4. Under "User variables", click "New"
  5. Variable name: `PATH`
  6. Variable value: `C:\Program Files\NCBI\blast\bin`
  7. Click OK and restart PowerShell

### "Docker not found" or "Docker daemon not running"
- Make sure Docker Desktop is installed and running
- Check: `docker --version` in PowerShell
- If not running, launch Docker Desktop from the Start menu

### "SibeliaZ failed" (when using Docker)
- Verify the image exists: `docker images | findstr sibeliaz`
- Try pulling/building again
- Check Docker Desktop is running and has sufficient resources (Settings → Resources)

### "permission denied" or "insufficient_scope" errors
- Make sure you're logged into Docker if using the pre-built image:
  ```powershell
  docker login
  ```
- Then pull the image again

## File Versions

The recommended versions above (BLAST+ 2.17.0) are current as of the time this guide was written. If newer versions are available at the NCBI FTP site, you can use them instead - just download the Windows AMD64 version (`.exe` or `-win64.tar.gz`).

## Getting Help

If you encounter issues:
1. Check that all dependencies are installed: `blastn -version` and `docker --version`
2. Verify your `config.txt` settings
3. Check the program's debug output for error messages
