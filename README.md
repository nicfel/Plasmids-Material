Github repository for the manuscript: 
# Quantifying plasmid movement in drug-resistant *Shigella* species using phylodynamic inference

**Authors**  
- Nicola F. Müller<sup>1,2,†</sup>  
- Ryan R. Wick<sup>3</sup>  
- Louise M. Judd<sup>4</sup>  
- Deborah A. Williamson<sup>5</sup>  
- Trevor Bedford<sup>2,6</sup>  
- Benjamin P. Howden<sup>3,4</sup>  
- Sebastián Duchêne<sup>3,7,‡</sup>  
- Danielle J. Ingle<sup>3,‡,†</sup>  

<sup>1</sup>Division of HIV, ID and Global Medicine, University of California San Francisco, CA, USA  
<sup>2</sup>Vaccine and Infectious Disease Division, Fred Hutchinson Cancer Center, Seattle, WA, USA  
<sup>3</sup>Department of Microbiology and Immunology at the Peter Doherty Institute for Infection and Immunity, The University of Melbourne, Melbourne, VIC, AUS  
<sup>4</sup>Centre for Pathogen Genomics, University of Melbourne (Doherty Institute)  
<sup>5</sup>Department of Infectious Diseases at the Peter Doherty Institute for Infection and Immunity, The University of Melbourne, Melbourne, VIC, AUS  
<sup>6</sup>Howard Hughes Medical Institute, Seattle, WA, USA  
<sup>7</sup>EDID unit, Department of Computational Biology, Institut Pasteur, Paris, France  

<sup>†</sup>Corresponding authors: nicola.felix.mueller@gmail.com, danielle.ingle@unimelb.edu.au  
<sup>‡</sup>These authors contributed equally to this work

---

## Abstract
The ‘silent pandemic’ of antimicrobial resistance (AMR) represents a significant global public health threat. These AMR genes are often carried on mobile elements, such as plasmids. The horizontal movement of these plasmids allows AMR genes, and the resistance they confer to key therapeutics, to disseminate throughout a population. However, quantifying the movement of plasmids remains challenging with existing computational approaches.

Here, we introduce a novel method for reconstructing and quantifying the movement of plasmids in bacterial populations over time. To do this, we model the co-evolution of chromosomal and plasmid DNA using a joint coalescent and plasmid transfer process in a Bayesian phylogenetic network framework. Our approach captures differences in the evolutionary histories of plasmids and chromosomes to identify instances where plasmids likely moved between bacterial lineages, while rigorously accounting for uncertainty in the data.

We apply this new method to a five-year dataset of *Shigella*, exploring transfer rates for five different plasmids that harbor various AMR and virulence profiles. In doing so, we reconstruct the co-evolution of the large *Shigella* virulence plasmid with the chromosomal DNA, quantify the higher transfer rates of three small plasmids that circulate among lineages of *Shigella sonnei*, and describe the recent dissemination of a multidrug-resistant plasmid between *S. sonnei* and *S. flexneri* lineages. Notably, this plasmid appears to have transferred in multiple independent events alongside an overall increase in its prevalence since 2010.

Our approach provides a powerful framework for understanding the evolutionary dynamics of plasmids carrying AMR genes as they are introduced, circulate, and are maintained in bacterial populations.
