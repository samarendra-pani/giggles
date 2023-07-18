# Giggles

Giggles is a software for genome-wide genotyping using long read sequencing data and a pangenome referenece.

## Contents

- [Installation](#installation)
    - [Building from Source using Pip](#building-from-source-using-pip)
- [Required Input Files](#required-input-files)
    - [Pangenome VCF](#pangenome-vcf)
        - [Creating a Pangenome VCF with assemblies](#creating-a-pangenome-vcf-using-assemblies)
    - [Alignments](#alignments)
    - [Haplotags](#haplotags)
- [Command Line Interface](#command-line-interface)
    - [Positional Arguments](#positional-arguments)
    - [Optional Arguments](#optional-arguments)
    - [Input pre-processing, selection and filtering](#input-pre-processing-selection-and-filtering)
    - [Realignment](#realignment)
    - [Emission Probability Parameters](#emission-probability-parameters)
    - [HMM Parameters](#hmm-parameters)
- [Test Examples](#test-example)


## Installation

### Building from source using Pip

```sh
git clone https://github.com/samarendra-pani/giggles.git
cd giggles
pip install -e .
```
This installs all the dependencies.

## Required Input Files

### Pangenome VCF

#### Creating a pangenome vcf using assemblies

### Alignments

### Haplotags


## Command Line Interface

```sh
usage: giggles genotype [-h] 
                        [-o OUTPUT] 
                        [--reference-fasta FASTA] 
                        [--rgfa rGFA] 
                        [--read-fasta FASTA] 
                        [--haplotag-tsv HAPLOTAG] 
                        [--max-coverage MAX_COV] 
                        [--mapping-quality QUAL] 
                        [--sample SAMPLE] 
                        [--chromosome CHROMOSOME] 
                        [--gt-qual-threshold GT_QUAL_THRESHOLD] 
                        [--keep-untagged] 
                        [--realign-mode MODE] 
                        [--overhang OVERHANG]
                        [--gap-start GAPSTART] 
                        [--gap-extend GAPEXTEND] 
                        [--mismatch MISMATCH] 
                        [--match-probability MATCH_PROBABILITY] 
                        [--mismatch-probability MISMATCH_PROBABILITY] 
                        [--insertion-probability INSERTION_PROBABILITY] 
                        [--deletion-probability DELETION_PROBABILITY] 
                        [--recombrate RECOMBRATE] 
                        [--eff-pop-size EFFPOPSIZE]
                        VCF 
                        [READS ...]
```
### Positional Arguments

- `VCF`: The pangenome VCF contains phased variants. This will be used as the reference panels. The output VCF will contain the same set of variants defined here. (refer [here](#pangenome-vcf))

- `READS`: BAM/GAF files containing the alignment information. (refer [here](#alignments))

### Optional Arguments

- `--reference-fasta`: The reference FASTA file has to be provided with BAM inputs.

- `--rgfa`: The reference GFA file has to be provided with GAF inputs. It should have BO and NO tag information (which can be done using [gaftools](https://github.com/marschall-lab/gaftools)).

- `--read-fasta`: This refers to the FASTA file which was aligned to the rGFA. This has to be provided with the GAF input (since GAF alignment lacks the alignment sequence).

- `--haplotag-tsv`: T

### Input Pre-processing, Selection and Filtering

- `--max-coverage`:

- `--mapping-quality`:

- `--sample`:

- `--chromosome`:

- `--gt-qual-threshold`:

- `--keep-untagged`:

### Realignment

- `--realign-mode`:

- `--overhang`:

- `--gap-start`:

- `--gap-extend`:

- `--mismatch`: 

### Emission Probability Parameters

- `--match-probability`: 

- `--mismatch-probability`: 

- `--insertion-probability`: 

- `--deletion-probability`: 

### HMM Parameters

- `--recombrate`: 

- `--eff-pop-size`: 

## Test Example





\pagebreak

# Developing

## Contents

- [Installation](#installation-1)
    - [Using Conda Environment](#conda-environment-installation-recommended-installation)
    - [Using Python Virtual Environment](#python-virtual-environment-installation)
- [Working with Cython Interface](#working-with-the-cython-interface)

## Installation

### Conda Environment Installation (Recommended Installation)

Giggles can be installed in a developmental environment using the following commands:

```sh
git clone https://github.com/samarendra-pani/giggles.git
conda create -n giggles-dev python=3.10
conda activate giggles-dev
cd giggles
pip install -e .
```

To use Giggles, user needs to activate the conda environment.
```sh
conda activate giggles-dev
```

To exit the conda environment, you can use `conda deactivate` or `conda activate base`.

Remove the conda environment to uninstall the local development version.

### Python Virtual Environment Installation

Giggles can be installed in a developmental environment using the following commands:

```sh
git clone https://github.com/samarendra-pani/giggles.git
cd giggles
python -m venv venv
source venv/bin/activate
pip install -e .
```

To use Giggles, user needs to activate the python virtual environment environment.
```sh
cd giggles
source venv/bin/activate
```

To exit the python virtual environment, you can use `deactivate`.

Remove the folder `venv` to remove the local development version

## Working with the Cython Interface

The Giggles genotyping algorithm is written in C++, as are many of the core data structures such as the “Read” class. To make the C++ classes usable from Python, we use Cython to wrap the classes. All these definitions are spread across multiple files. To add new attributes or methods to an existing class or to add a new class, changes need to be made in different places.

Let us look at the “Read” class. The following places in the code may need to be changed if the Read class is changed or extended:

- ``src/read.cpp``: 
  Implementation of the class (C++).

- ``src/read.h``: 
  Header with the class declaration (also normal C++).

- ``giggles/cpp.pxd``: 
  Cython declarations of the class. This repeats – using the Cython syntax this time – a subset of the information from the ``src/read.h`` file. This duplication is required because Cython cannot read ``.h`` files (it would need a full C++ parser for that).

  Note that the ``cpp.pxd`` file contains definitions for *all* the ``.h``
  headers. (It would be cleaner to have them in separate ``.pxd`` files, but
  this leads to problems when linking the compiled files.)

- ``giggles/core.pxd``: 
  This contains declarations of all *Cython* classes wrapping C++ classes. Note that the class ``Read`` in this file has the same name as the C++ class, but that it is not the same as the C++ one! The distinction is made by prefixing the C++ class with ``cpp.``, which is the name of the module in which it is declared in (that is, the C++ class ``Read`` is declared in ``cpp.pxd``). The wrapping (Cython) class ``Read`` stores the C++ class in an attribute named ``thisptr``. If you add a new class, it needs to be added to this file. If you only modify an existing one, you probably do not need to change this file.

- ``giggles/core.pyx``: 
  The Cython implementation of the wrapper classes. Again, the name ``Read`` by itself is the Python wrapper class and
  ``cpp.Read`` is the name for the C++ class.
  
  Before adding yet more C++ code, which then requires extra code for wrapping it, consider writing an implementation in Cython instead. See ``readselect.pyx``, for example, which started out as a Python module and was then transferred to Cython to make it faster. Here, the Cython code is not merely a wrapper, but contains the implementation itself.