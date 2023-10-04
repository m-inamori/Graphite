# Graphite

GRAPHITE
GRAph used tool to correct, PHase, and ImpuTE VCF


## Usage

graphite VCF ped [-m map] [-t num_threads] [-c chrom indices] [-p lower progenies] out.
chrom indices: (index|first:last)[,(index|first:last),..]

*VCF*             : input VCF

*ped*             : pedigree file
                    A text file with 4 columns and without header.
                    Its delimiter is a space.
                    Its columns are family name, sample name, a parent name, and other parent name.
                    Parent order is reflected in the VCF's Genotype.

*out*             : output VCF

### option

*map*             : genetic map file
                    A CSV file with 3 columns and without header.
                    Its deliimiter is a comma.
                    Its columns are scaffold name, cM, and Mbp.

*chrom indices*   : output scaffolds' indices

*lower progenies* : the lower limit on the number of progeny that can be regarded as large families.

#### 

## License
