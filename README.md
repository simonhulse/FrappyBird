<!--README.md-->
<!--Simon Hulse-->
<!--simonhulse@protonmail.com-->
<!--Last Edited: Wed 26 Feb 2025 05:03:12 PM EST-->

# FrappyBird

Script for processing FRAP data.

## Usage

### Filename requirement

**TL;DR:** *To make this script work, make sure the `.csv` files you wish to
process follow the same naming convention as the files in this repository. For
example:* `r4-850uM_pu30-673M_50min A - 2_.csv`

For this script to work, it is necessary for all the `.csv` files that you wish
to process to match the following regular expression:

    .* [A-Z].*\.csv

This means it should consist of:

1. Anything you want (i.e. experiment conditions), followed by a space
2. A capital letter after the space
3. Anything you want after the capital letter
4. A `.csv` extension

### Running the processing script

Load `MATLAB`, and enter the following commands at the prompt:

    >> cd DIRECTORY/CONTAINING/THE/CSV/FILES/
    >> !curl -o FrappyBird.m https://raw.githubusercontent.com/simonhulse/FrappyBird/refs/heads/main/FrappyBird.m
    >> FrappyBird

The result of the processing is stored in the newly-formed `output` directory.
