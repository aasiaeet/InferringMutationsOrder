---
title: Inferring Mutation Order (IMO)
author: Kevin R. Coombes
date: 2018-05-02
---

The ultimate goal of this project is to develop better methods to
infer the most likely order of mutations in different kinds of
cancer. 

The initial idea is that existing data sets describing the
co-occurrence patterns of mutations can be used to infer mutation
order.  In the simplest case, if the only time we see mutation B in a
patient's tumor is when we also see mutation A, but we see some cases
of mutation A without B, then it is pretty asy to conclude that A
occurred before B. In the real world, things are always more
complicated, and various models have been trying to extend this simple
idea to something that might apply to real data.

# Project 1: Literature Review

1. We will start by reviewing the existing methods.
2. We will collect existing code or, if necessary, possibly
   reimplement the method.
3. We wil apply all methods to the same data sets:
    a. TCGA colon (COAD) and rectal (READ) cancers
	b. TCGA cutaneous (skin) melanoma (SKCM)
4. We will then write up the results, reporting on which methods seem
   to work better.

Note that the chosen cancer types are characterized by the fact that,
because of extensive screening programs, many precancerous lesions
have alrady been studied independently of TCGA. As a result, there are
existing models of carcinogenesis that are based on actual data.

# Project 2: Incorporate MAF

All of the existing methods treat mutations as a binary
(absent/present) event. But sequencing studies always yield a
"mutation allele frequency" (MAF) that describes the percentage of
cells that have the mutation. We hope to incorporate the MAF intot he
analysis, based on the "prior" belief that mutations that occur in
more cells usully occur earlier than mutations that occur in fewer
cells.

# Repository Structure

We structure the github repository so it contains three primary
components: 

1. `code` contains all the computer code to perform the analysis
2. `doc` contains both the documentation for how the code works as
   well as any manuscripts that 
3. `results` contains any results of the anlaysis that should be saved
   as support for resiulting papers and presentations. 

In addition, there are three structured components that are part of
the project but not maintained under version control. We seaprate
these portions because that are usually large, static, and binary,
making them unsuitable for a text-based code repository.  The parts
are 

1. `raw` data: this folder contains the primary source files
   containing the data we use for the analysis.
2. `clean` data: this folder contains any processed data used as part
   of the analysis. In moset cases, the `clean` data is produced by
   doing something to the `raw` data (such as normalization), and this
   processing is docuemnted in the earliest scripts in the `code`
   folder.
3. `scratch`: this folder contains intermediate results. They can be
   regenerated by the existing `code` but, because the process is
   often time-consuming, we prefer to save the intemediate results
   somewhere. The distinction between `scratch` and `results` is
   whether or not the items inside them are intended to be used in
   final papers or presentations.

## Indirection via JSON
We arrange things so that individuals can store their data and scratch
space wherever they want. (On their local hard drive, a network drive,
somehwere in the cloud; dropbox, etc.) But, we want the code to run
without editing on everybnody's individual machihne. So, we adopt the
following convention. Each user must create a file called
	`$HOME/Paths/imo.json`
Here, as usual, `$HOME` refers to the user's home directory on the
local machine. The subfolder `Paths` is hard-coded and cannot be
changes. The file name `imo.json` is specific to this project; other
projects will use different JSON files stored inthe same directory.

A sampe `imo.json` file is provided at the top-levle of the
repository. Yu can copy this into the correct place and edit it to
point to the desired locations on your machine.
