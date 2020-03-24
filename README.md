# Illumina-CoVID19-RefMap

### Remove PCR Primers

Illumina sequencing is done using PCR amplicons. We have to remove PCR primers before doing any processing of the data. Depending on the primer scheme,
either we have [V1](https://raw.githubusercontent.com/artic-network/artic-ncov2019/master/primer_schemes/nCoV-2019/V1/nCoV-2019.tsv)  primers or [V2](https://raw.githubusercontent.com/artic-network/artic-ncov2019/master/primer_schemes/nCoV-2019/V2/nCoV-2019.tsv) primers in the fastq sequences.

The program cutadapt (with -b option for primers) is used for removing primers. There is an automated script available on the alpha server to process fastq files. It can be used for V1 and V2 primer schemes. If a new primer scheme is used, those primer sequences should be added to this script. 

Usage is:
```
/home4/nCov/Sreenu/Scripts/removePrimers.sh file1.fastq file2.fastq
```

This will remove the primers and store the output in "file1-woP.fq" and "file2-woP.fq". Total number of reads before and after primer timming should be same.

### Quality trimming 

Then we are trimming the low quality reads and retaining reads with more than 75bp using trim_galore program.
```
trim_galore -length 75 -q 30 --illumina --paired file1-woP.fq file2-woP.fq
```

Results will be stored in "file1-woP_val_1.fq" and "file2-woP_val_2.fq"

### Mapping

Trimmed reads are mapped using tanoti to Wuhan nCov genome (MN908947.3). This file is stored in "/home4/nCov/Sreenu/Ref/MN908947.3.fa"

```
tanoti -p 1 -r /home4/nCov/Sreenu/Ref/MN908947.3.fa -i file1-woP_val_1.fq file2-woP_val_2.fq -o file.sam 
```

### Mapping statistics and consensus sequence generation

Consensus sequence from file.sam is generate using SAM2CONSENSUS with default values. By default it generates a majority consensus. For true consensus sequence use "-c 50". It mimics reference genome and ignores INDELS.
Indel option will be added in next version. 

```
SAM2CONSENSUS -i file.sam -o file-con.fa
```

Assembly statistics are generated using  SAM_STATS and coverage is plotted using SameerReport 
```
SAM_STATS_tbl file.sam > file.stats 
SameerReport-Small file.sam
```

### Data for Climb upload

For phylogenetic analysis, ARTIC consortium is asking consensus sequence and only mapped reads in a sorted bam file. 

We have to strip LabIDs from file names and also from sam files. In below command we are removing comment line from sam header i.e "@CO:", printing only mapped reads and converting them to sorted bam file.

Sam fie comment header has command line arguements which sometimes contain original fastq file name with LabIDs. 

```
grep -v "^@CO" file.sam|samtools view -F4 -bS  - |samtools sort  -o  file.bam 
```

Finally, prepare the data for  CLIMB upload. LabIDs from the file names should be stripped before uploading. Consensus sequence and mapped reads in sorted bam format should be stored in respective directories.
I have created a script to copy and process data for upload. Pease change the code as per your directory structure

In below example. My analysis directory is "/home4/nCov/Sreenu/batch7b" and CLIMB upload directory is "/home4/nCov/Sreenu/ClimbData/CVR-Illumina-Batch7b/"


```
cd /home4/nCov/Sreenu/batch7b
find . -name "*.bam" -exec cp {} /home4/nCov/Sreenu/ClimbData/CVR-Illumina-Batch7b/. \;
find . -name "*-con.fa" -exec cp {} /home4/nCov/Sreenu/ClimbData/CVR-Illumina-Batch7b/. \;
cd  /home4/nCov/Sreenu/ClimbData/CVR-Illumina-Batch7b/
cp /home4/nCov/Sreenu/Scripts/renameConsensus.sh .
./renameConsensus.sh
rm renameConsensus.sh
```

This copies bam and consensus files to  /home4/nCov/Sreenu/ClimbData/CVR-Illumina-Batch7b folder, renames them and chnages the consensus fasta header appropriately

Upload this entire folder to CLIMB. [Instruction can be found here]((https://docs.covid19.climb.ac.uk/upload-instructions)
