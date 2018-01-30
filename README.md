# mgpfusion
mGPfusion is a Gaussian process based method for predicting stability changes upon single and multiple mutations  of proteins that complements the available experimental data with large amounts of simulated data. 
# mGPfusion

mGPfusion is a Gaussian process based method for predicting stability changes upon single and multiple mutations 
of proteins that complements the available experimental data with large amounts of simulated data. 
Our Bayesian data fusion model re-calibrates the experimental and in silico data sources and then learns 
a predictive model from this combined data. Our model is protein-specific and requires experimental data only 
regarding the protein of interest. mGPfusion has been developed at Aalto University.
* For a comprehensive description of mGPfusion see \[1\] 
* For an example usage of mGPfusion, see example.m

## Data
Our example data contains information for 15 proteins: pdb-structures from [The Protein data bank](http://rcsb.com), experimentally measured ddG-values 
from [Protherm](http://www.abren.net/protherm/) \[2\], and ddG-values simulated with [Rosetta](https://www.rosettacommons.org/) \[3\].

### ddg-values in csv-files and cell arrays (mat-files)
The first column in the cell arrays (second column in csv-files) contains the mutations defined using 
consecutive numbering. The third column in csv-files uses Rosetta's numbering. The signs of the ddG-values 
simulated with Rosetta have been changed to match the notation used for the muations from Protherm.
The scaling suggested by \[4\], `0.57 * y_S`, can be applied if the values are used as they are.
ddg_protherm contains experimentally measured ddG-values gathered from Protherm, ddg_rosetta_single contains
ddG-values for all single mutations simulated with Rosetta, and ddg_rosetta_multi contains ddG-values for all
multiple mutations simulated with Rosetta.

#### ddg_protherm
The experimentally measured ddG-values have been obtained from Protherm. If multiple measurements for the ddg-change upon a mutation were provided, an average was taken. 

#### Mutation numbering
We use consecutive numbering for mutations, which means that the first residue in the amino acid sequence 
has number 1, second 2, and so on, regardless if the residues are present in the pdb structure. For our 15
example proteins, this corresponds to the numbering used in the pdb structures, but this may not be the case 
for all proteins.

The numbering that is used by Rosetta gives numbers only to residues that are present in the pdb-structure.
This numbering often differes from our consecutive numbering and that used with pdb-structures.

## Copyright

### ddg_protherm.csv and ddg_protherm.mat
This data is licensed under [Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License](http://creativecommons.org/licenses/by-nc-sa/3.0/)

### minFunc
For the minFunc software, see the copyright.txt file in repository code/minConf

### Everything else in this mgpfusion-repository
Copyright 2017-2018 Emmi Jokinen

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use these files except in compliance with the License.
   You may obtain a copy of the License at
   
      http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

## References
\[1\] Jokinen, E., Heinonen, M., and Lähdesmäki, H. *mGPfusion: Predicting protein stability changes with Gaussian process kernel learning and data fusion*. (submitted)

\[2\] Kumar MD, Bava KA, Gromiha MM, Parabakaran P, Kitajima K, Uedaira H, Sarai A. (2006). [*ProTherm and ProNIT: thermodynamic databases for proteins and protein-nucleic acid interactions*](http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=16381846&dopt=Abstract). Nuleic Acids Res., 34:D204-6, Database issue

\[3\] Leaver-Fay, A., Tyka, M., Lewis, S. M., Lange, O. F., Thompson, J., Jacak, R.,
Kaufman, K., Renfrew, P. D., Smith, C. A., Sheffler, W., et al. (2011). *ROSETTA3:
an object-oriented software suite for the simulation and design of macromolecules*.
Methods in enzymology, 487, 545

\[4\] Kellogg, E. H., Leaver-Fay, A., and Baker, D. (2011). *Role of conformational sampling in computing mutation-induced changes in protein structure and stability*. Proteins: Structure, Function, and Bioinformatics, 79(3), 830-838.

