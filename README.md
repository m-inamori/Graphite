# Graphite

GRAPHITE
GRAph used tool to correct, PHase, and ImpuTE VCF

## Usage

  graphite [options]

  -i &lt;path&gt;/--input &lt;path&gt;     Input VCF file  
  -p &lt;path&gt;/--pedigree &lt;path&gt;  Pedigree file  
  -o &lt;path&gt;/--output &lt;path&gt;    Output file  

### Options

  -m &lt;path&gt;/--map &lt;path&gt;       Input map file (CSV format with 3 columns: scaffold name, cM, Mbp)  
  -r &lt;path&gt;/--ref &lt;path&gt;       Reference VCF file  
  -t &lt;int&gt;/--num-threads &lt;int&gt; Number of threads (default: 1)  
  --lower-progs &lt;int&gt;          Lower limit on the number of progenies considered as a large family  
  --not-impute-isolated        Do not impute isolated samples  
  --not-correct-isolated       Do not correct wrong genotypes of isolated samples  
  --out-isolated               Output isolated samples separately (valid only with --not-impute-isolated)  
  -c &lt;int&gt;/--chrom &lt;int&gt;       Impute only the chromosome specified by the 0-based index  
  --precision-ratio &lt;float&gt;    Control runtime for small pedigree HMM analysis (default: 1.0; larger values increase runtime)  
  --fast                       Shortcut for --precision-ratio=0.1 (optimized for faster runtime, reduced precision)  
  --precision                  Shortcut for --precision-ratio=10.0 (enhanced precision, increased runtime)  

#### 

## License
MIT License

Copyright (c) 2023 Minoru Inamori

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
