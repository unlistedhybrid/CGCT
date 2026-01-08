# CGCT: Comparative Genomics Circular Toolkit

A desktop application for comparative whole-genome analysis with interactive circular visualization and quantitative metrics.

## Overview

CGCT integrates whole-genome alignment with fine-scale homology detection and BLAST homology searching to provide a unified view of genomic conservation and divergence. The application combines **progressiveMauve** (large-scale synteny), **SibeliaZ** (local homology), and **BLAST+** (gene homology) and presents results in an interactive circular map with publication-quality output.

Built in Python 3.12 with Tkinter and Matplotlib, CGCT runs natively on Windows, macOS, and Linux.

## Key Features

### Hybrid Alignment Pipeline
- **progressiveMauve 2.5.0**: Identifies Locally Collinear Blocks (LCBs) free from rearrangement with gapped refinement
- **SibeliaZ integration**: Detects local homology (≈50 bp) in diverged or non-coding regions
- **Automatic merging**: Combines both outputs to capture structural context and fine-scale divergence
- **Deterministic alignment**: Produces identical results across all platforms

### Quantitative Analysis
- **Identity profiles**: Sliding-window nucleotide identity with coding/non-coding breakdown
- **ANIb calculation**: Fragment-based average nucleotide identity with coverage metrics
- **dDDH estimation**: Digital DNA-DNA hybridization from alignment statistics
- **GC analysis**: GC content and skew computed in real-time with adjustable window size

### Interactive Visualization
- **Circular maps**: Concentric tracks for genes, BLAST hits, identity, and GC metrics
- **Navigation tools**: Pan, zoom, rotation, box-zoom, and region-of-interest highlighting
- **Multi-format export**: PNG, SVG, PDF, TIFF, JPG, or EPS up to 19200×10800 pixels
- **Real-time updates**: Adjust visualization parameters without re-running alignments

### Flexible Input
- **Genomic sequences**: FASTA format (single or multi-contig)
- **Annotations**: GFF or GenBank format with automatic coordinate conversion
- **Regions of interest**: BED format for highlighting specific loci
- **Pre-computed alignments**: Supply Mauve (XMFA) or SibeliaZ (MAF) files directly
- **Automatic validation**: Cross-check sequence and annotation boundaries before analysis

### Comprehensive Output
- **Publication-quality figures**: Circular maps in multiple formats at user-defined resolution
- **Excel workbook**: Multi-sheet statistics including all alignment blocks, variants, and metrics
- **FASTA extraction**: Export any genomic region for downstream PCR or synthesis

## Installation

Download the appropriate version for your platform from [Releases](https://github.com/unlistedhybrid/CGCT/releases):

- **Windows**: `CGCT-windows-x64.zip`
- **Linux**: `CGCT-linux-x64.tar.gz`
- **macOS Intel**: `CGCT-macos-x64.tar.gz`
- **macOS ARM64**: `CGCT-macos-arm64.tar.gz`

Extract and run directly—no installation required. Each platform bundle includes:
- Pre-compiled progressiveMauve 2.5.0
- Platform-specific README with installation and troubleshooting information

### System Requirements

- **Windows**: Windows 10 or later
- **macOS**: macOS 10.13 or later
- **Linux**: Ubuntu 20.04 or later (or equivalent)

### SibeliaZ on Windows and macOS

CGCT can run SibeliaZ on Windows and macOS through containerization:
- **Windows**: WSL2 is recommended. A `.bat` installer script is included to automate WSL2 command-line setup. Alternatively, Docker can be used with the provided Dockerfile or by pulling a pre-built image: `docker pull aaronmauve/sibeliaz-conda:latest`
- **macOS**: A Dockerfile is provided, or pull the pre-built image: `docker pull aaronmauve/sibeliaz-conda:latest`

See the platform-specific README for detailed setup instructions.

## Usage

1. **Open CGCT** and load a query genome (FASTA) and reference genome (FASTA)
2. **Add annotations** (optional): Supply GFF or GenBank files for coding region analysis
3. **Run alignment**: Select "Mauve only", "SibeliaZ only", or "Hybrid" mode
4. **Explore results**: Navigate the circular map, adjust visualization parameters, and hover over regions for details
5. **Export**: Save figures in your preferred format and statistics as Excel workbook

For platform-specific setup and troubleshooting, see the README included with your download.

## Technical Details

### Architecture
- **Language**: Python 3.12
- **GUI**: Tkinter
- **Graphics**: Matplotlib, DNA Features Viewer
- **Genomics**: Biopython, NumPy, Pandas
- **Alignment**: progressiveMauve 2.5.0, SibeliaZ

### Dependencies

CGCT calls progressiveMauve and SibeliaZ as external executables. Pre-compiled progressiveMauve 2.5.0 is included in the distribution.

## License

This project is licensed under the **Business Source License 1.1 (BSL)**, with automatic conversion to **GPLv3 on January 7, 2029**.

### Summary
- **Non-production use**: Free (development, testing, academic research)
- **Production use**: Requires a commercial license
- **BSL 1.1 terms**: See [LICENSE.txt](LICENSE.txt)

### Dependencies
- **progressiveMauve**: GNU General Public License v2 (see [COPYING](COPYING))
- **Python packages**: MIT, BSD, PSF, and Biopython licenses (see [LICENSE.txt](LICENSE.txt))

For commercial licensing or production use inquiries, contact: `244086359+unlistedhybrid@users.noreply.github.com`

## Citation

If you use CGCT in published research, please cite:

```
Jackson, A., Thomas, M. (2026). CGCT: Comparative Genomics Circular Toolkit. 
https://github.com/unlistedhybrid/CGCT
```

## Authors

**Aaron Jackson** ([@unlistedhybrid](https://github.com/unlistedhybrid))  
**Michael Thomas** ([ORCID: 0000-0003-4646-5324](https://orcid.org/0000-0003-4646-5324))

## Contributing

We welcome bug reports, feature requests, and contributions. Please open an [Issue](https://github.com/unlistedhybrid/CGCT/issues) or submit a Pull Request.

## Support

- **Documentation**: See platform-specific README in your bundle
- **Issues**: [GitHub Issues](https://github.com/unlistedhybrid/CGCT/issues)
- **progressiveMauve**: [progressiveMauve Repository](https://github.com/unlistedhybrid/progressiveMauve-modernization)
- **SibeliaZ**: [SibeliaZ Documentation](https://github.com/medvedevgroup/SibeliaZ)
