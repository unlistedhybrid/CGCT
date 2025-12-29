# CGCT Installation Guide - Linux AMD64

This guide covers installing CGCT on Linux (AMD64 architecture).

## Prerequisites

- Linux with AMD64 processor (x86_64)
- Root/sudo access for installing software
- Common distributions: Ubuntu, Debian, CentOS, Fedora, etc.

## Step 1: Install BLAST+

BLAST+ is required for sequence alignment and homology detection.

### Download BLAST+

Visit: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

Choose based on your package manager:

**For Debian/Ubuntu/Linux Mint:**
- Download: **`ncbi-blast-2.17.0+-x64-linux.tar.gz`** (Recommended)

**For CentOS/RHEL/Fedora (RPM-based):**
- Download: **`ncbi-blast-2.17.0+-2.x86_64.rpm`** (Recommended)

### Install BLAST+ (tar.gz - Debian/Ubuntu/Others)

1. Download the `.tar.gz` file
2. Open Terminal
3. Navigate to your Downloads folder:
   ```bash
   cd ~/Downloads
   ```
4. Extract the archive:
   ```bash
   tar -xzf ncbi-blast-2.17.0+-x64-linux.tar.gz
   ```
5. Move it to a standard location:
   ```bash
   sudo mv ncbi-blast-2.17.0+ /usr/local/ncbi-blast
   ```
6. Add it to your PATH:
   ```bash
   echo 'export PATH="/usr/local/ncbi-blast/bin:$PATH"' >> ~/.bashrc
   source ~/.bashrc
   ```
   Or if using zsh:
   ```bash
   echo 'export PATH="/usr/local/ncbi-blast/bin:$PATH"' >> ~/.zshrc
   source ~/.zshrc
   ```

### Install BLAST+ (RPM - CentOS/RHEL/Fedora)

1. Download the `.rpm` file
2. Install it:
   ```bash
   sudo rpm -i ncbi-blast-2.17.0+-2.x86_64.rpm
   ```

### Verify Installation

Open Terminal and run:
```bash
blastn -version
```

You should see output showing the BLAST version.

## Step 2: Install SibeliaZ Natively

SibeliaZ can be installed directly on Linux without any virtualization or containers.

### Install Conda (if not already installed)

1. Download Miniconda:
   ```bash
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   ```

2. Run the installer:
   ```bash
   bash Miniconda3-latest-Linux-x86_64.sh
   ```

3. Follow the installer prompts and restart your terminal

### Install SibeliaZ via Conda

1. Create a conda environment:
   ```bash
   conda create -y -n sibeliaz python=3.11
   ```

2. Activate the environment:
   ```bash
   conda activate sibeliaz
   ```

3. Install SibeliaZ and twopaco:
   ```bash
   conda install -y -c bioconda sibeliaz twopaco
   ```

4. Test the installation:
   ```bash
   sibeliaz --version
   ```

   You should see the SibeliaZ version printed.

## Step 3: Configure CGCT

Edit the `config.txt` file (in the same directory as the CGCT program).

The default settings should work fine for Linux with native SibeliaZ installation.

## Running CGCT

1. Activate the SibeliaZ conda environment before launching CGCT:
   ```bash
   conda activate sibeliaz
   ```

2. Launch CGCT from the same terminal

3. The program will automatically find:
   - BLAST+ in `/usr/local/ncbi-blast/bin`
   - SibeliaZ in the activated conda environment

## Troubleshooting

### "blastn not found"
- Check installation: `blastn -version`
- Verify PATH:
  ```bash
  echo $PATH
  ```
- Should include `/usr/local/ncbi-blast/bin`
- If missing, re-run the PATH export command and restart your terminal

### "sibeliaz not found"
- Make sure conda environment is activated:
  ```bash
  conda activate sibeliaz
  ```
- Verify installation:
  ```bash
  sibeliaz --version
  ```
- If not found, reinstall SibeliaZ:
  ```bash
  conda install -y -c bioconda sibeliaz twopaco
  ```

### Conda activation doesn't work
- Verify conda is installed:
  ```bash
  conda --version
  ```
- If conda command not found, check the installer output - you may need to close and reopen your terminal
- You can also initialize conda manually:
  ```bash
  ~/miniconda3/bin/conda init
  ```

### "python" or other dependency errors
- Make sure to activate the sibeliaz environment:
  ```bash
  conda activate sibeliaz
  ```
- The environment includes Python 3.11 and all required dependencies

## File Versions

The recommended versions above (BLAST+ 2.17.0) are current as of this guide. Check the NCBI FTP site for newer versions and download the `x64-linux.tar.gz` version (or RPM if using an RPM-based distribution).

## Performance Notes

- Native SibeliaZ installation on Linux is very fast
- No virtualization or emulation overhead
- All tools run natively on your system

## Getting Help

If you encounter issues:
1. Verify BLAST+: `blastn -version`
2. Verify SibeliaZ: 
   ```bash
   conda activate sibeliaz
   sibeliaz --version
   ```
3. Check `config.txt` is in the correct location
4. Ensure the SibeliaZ conda environment is activated when launching CGCT
5. Review the program's error messages for specific details

