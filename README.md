SAIGE - FAMILY :

Software depndencies: 
Nextflow (https://www.nextflow.io/docs/latest/install.html)
Singularity (https://hub.docker.com/r/pennbiobank/saige)

Getting started (on LPC):
1.Make sure singularity and nextflow are loaded pipeline will not run otherwise
module load nextflow
module load singularity

2. Recreate SAIGE FAMILY using symbolic links for every directory EXCEPT the configs this should be copied
   SAIGE_FAMILY/
       | configs/
       | workflows/
       | processes/
       | scripts/
       | nextflow.configs
3. In the config for your pipeline of your choice (EXWAS) edit the corresponding configuration file with desired file paths for analyses
4. Run using cluster profile
    nextflow run ./workflows/saige_exwas.nf -profile cluster
    nextflow run ./workflows/saige_exwas.nf -profile cluster -resume
   



