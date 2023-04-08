#SEEC-nt
Protein sequence evolution with Direct Coupling Analysis (DCA) in Matlab.

#Overview
Given an amino acid sequence and its associated nucleotide sequence, the package produces an evolutionary trajectory of sequences using constraints generated from the eij (coupings) and hi (local fields) parameters from the Potts model of the protein family. The user provides the coupling and local field inputs by inferring them using any DCA method of their choosing.  They can also control the temperature (T) of the conditional probability distribution sampling step, which impacts sequence divergence, and they can also specify the number of evolutionary steps in the simulation.  

#Features

This package is similar to our original model (https://github.com/AlbertodelaPaz/SEEC) but with some important changes to make the simulation more biologically relevant, as was previously suggested by Basardi et al (https://doi.org/10.1093/molbev/msab321) . Whereas previously, any sampled amino acid could replace the existing amino acid, here, the movement in sequence space is restricted to amino acids that are one base change away from the codon encoding the existing amino acid.    For example, say the chosen site has Ser with the associated codon being UCU; if Tyr is sampled, it would be accepted, since there is a codon (UAU) that differs by a single nucleotide.  If, on the other hand, the Ser had UCA as its current codon, then Tyr wouldn't be accepted since Tyr is only coded with UAU and UAC and both differ by two nucleotides.

In the case where UCU is coding for Ser and Ser is sampled again, then we would choose at random between UCU, UCC, UCA and UCG since all of them differ by 0 or one nucleotide, but we would ignore AGU and AGC since they are at 2 and 3 nucleotides of distance, respectively.  In this way, synonymous mutations can be made.

In our original SEEC model, the length of the input sequence could be lengthened or shortened during the simulation due to replacement of gaps with amino acids or the replacement of amino acids with gaps. Here, instead of being able to choose any site for potential substitution, the length of the protein is maintained throughout the simulation, as positions with existing gaps cannot be chosen for substitution, and if a a position with an amino acid is chosen and a gap is sampled from the conditional probability distribution, the distribution is re-sampled until an amino acid is chosen.  If a gap is sampled 100 times, the current amino acid and its codon remain, a new site is chosen and a new step is registered.  


#Install
To install, download the package into the main MATLAB folder.

#Technologies
Project is created with:
* Matlab version: 2020b 

Probevolution_locality.m main script to run the evolutionary simulation.
Randvar.m  returns a realization for a discrete random variable given its cumulative probability distribution.
completeDCA.m infers Potts model parameters using mean-field DCA.  Outputs parameters, the input alignment in number format, and a list of DI and MI values for each pair of positions i and j.
getAlignmentFromFile.m  converts amino acid and nucleotide .fasta files into number format.
Molecule2Num.m  convert .fasta files to number format.
Generalhamiltonian.m Computes the Potts Hamiltonian for each amino acid sequence in a provided array. 
siteprobdistribution.m  Calculates the conditional probability distribution for sites. Output is used for the Probevolution_locality script. 
isNB.m uses the standard genetic code to determine if the amino acid that was sampled from the distribution is a neighbor (one base away) of the current amino acid.  


#Usage
Open Matlab and enter the package folder as the working directory.  
Ensure the amino acid and nucleotide sequences (both .fasta format) are also in this folder.
Infer family Potts model parameters using CompleteDCA.m, which uses mfDCA, or use another method of DCA.  
Convert nucleotide sequence to number format using getAlignmentFromFile.m.
Run Probevolution_locality, inputting any positive scalar for T.  
function [Trajectory_amino,Trajectory_nucleo, H, sitecount,sustcount,Timeline,Sustitutiontime,R,Sustitutiontimeline,mutatedsites]=Probevolution_locality(S_amino,S_nucleo,e,h,T,M)


#Acknowledgements
Any publication resulting from applications of SEEC-nt should cite:
Sophia Alvarez*, Charisse M. Nartey*, Nicholas Mercado, Alberto de la Paz, Tea Huseinbegovic, Faruck Morcos.  “Evolutionary simulation under epistatic constraints yields trajectories with in vivo functionality.”  Under revision.

