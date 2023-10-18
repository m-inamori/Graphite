# Graphite

GRAPHITE
GRAph used tool to correct, PHase, and ImpuTE VCF


## Usage

graphite  -i VCF -p ped [-m map] [-t num_threads] [-f family indices] [-c chrom indices] [--lower-progs lower num progenies] [--large-only] [--not-impute-isolated [--out-isolated]] -o out.
family indices: (index|first:last)[,(index|first:last)[,..]]
chrom indices: same as family indices.

*VCF*             : input VCF

*ped*             : pedigree file
                    A text file with 4 columns and without header.
                    Its delimiter is a space.
                    Its columns are family name, sample name, a parent name, and other parent name.
                    Parent order is reflected in the VCF's Genotype.

*out*             : output VCF

### optional

*map*             : genetic map file
                    A CSV file with 3 columns and without header.
                    Its deliimiter is a comma.
                    Its columns are scaffold name, cM, and Mbp.
                    If not specified, 1cM = 1Mbp.

*chrom indices*   : output scaffolds' indices

*lower num progenies* : the lower limit on the number of progeny that can be regarded as large families.

--large-only      : output large families only

--not-impute-isolated : not impute isolated samples

--out-isolated    : output not imputed isolated samples

#### 

## License
