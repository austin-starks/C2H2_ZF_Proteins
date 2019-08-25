HOW TO USE:

To use this code, you must have the following requirements:
   Python 3+
   HMMER 2.3.2

Download the code. In the same folder that contains "main", create a folder with the
following names:
   analysis
   chromosomes
   hmmer
   output
   
Databases:
 There should be a folder named databases present. Unzip the folder.
 
HMMER:
 There should be a folder named HMMER present. Unzip the folder.
 
Analysis:
 This folder is used to compare the predicted change in DNA-binding specificity.
 To do this, visit zf.princeton.edu and click "Predict PWS". Type in the ZF protein
 sequence and click send. Click on the appropiate ZF domain and chaneg the prediction 
 model to be Polynomnial SVM. Click submit and download the PWM. Save the file in the 
 following format:
 
 GeneName_DomainNumber_o.txt - for an original gene.
 GeneName_DomainNumber_m.txt - for a mutated gene.
 
Chromosomes:
 This folder is used to store the chromosomes. Chromsomes can be downloaded from: 
 ftp://ftp.ensembl.org/pub/release-97/fasta/homo_sapiens/dna/
 

Finally, run project.py as a script to get the main output files. There are other 
functions in the file that can do a variety of things. Feel free to explore. If 
you have any questions, please email me at as2779[at]cornell.edu.
 
