###For the samples from NEB library
#unify name
for f in *_L001_R1_001.fastq.gz; do mv -- "$f" "${f%_L001_R1_001.fastq.gz}_R1.fastq.gz"; done
for f in *_L001_R2_001.fastq.gz; do mv -- "$f" "${f%_L001_R2_001.fastq.gz}_R2.fastq.gz"; done

#concat (put forward and reverse sequence together)
cd ../
mkdir all_concat
gunzip copy/*.gz

for i in copy/*R1.fastq
do
R1=$i;
R2=${R1%_R1.fastq}_R2.fastq;
SAMPLE=`basename ${R1%_R1.fastq}`;
echo $SAMPLE
echo $R1
echo $R2

cat $R1 $R2 > all_concat/$SAMPLE.catF.fastq
cat $R2 $R1 > all_concat/$SAMPLE.catR.fastq
done

#good orientation
mkdir good_orientation
for i in all_concat/*catF.fastq; 
do R1=$i; R2=${R1%catF.fastq}catR.fastq; 
SAMPLE=`basename ${R1%.catF.fastq}`; 
echo $SAMPLE 
cutadapt -g CCTACGGGNGGCWGCAG -G GACTACHVGGGTATCTAATCC -o good_orientation/$SAMPLE.orient_R1.fastq -p good_orientation/$SAMPLE.orient_R2.fastq $R1 $R2 --discard-untrimmed; 
done

gzip good_orientation/*.fastq

conda deactivate
exit

###For the samples from 2-step PCR
cd (file position)

conda activate /home/ke36dar/data/programs/cutadaptenv

# unify names
for f in *_L001_R1_001.fastq.gz; do mv -- "$f" "${f%_L001_R1_001.fastq.gz}_R1.fastq.gz"; done
for f in *_L001_R2_001.fastq.gz; do mv -- "$f" "${f%_L001_R2_001.fastq.gz}_R2.fastq.gz"; done

# trim
mkdir trimmed
for i in *_R1.fastq.gz; 
do R1=$i; R2=${R1%_R1.fastq.gz}_R2.fastq.gz; 
SAMPLE=`basename ${R1%_R1.fastq.gz}`; 
echo $SAMPLE 
cutadapt -g CCTACGGGNGGCWGCAG -G GACTACHVGGGTATCTAATCC -o trimmed/$SAMPLE.trim_R1.fastq.gz -p trimmed/$SAMPLE.trim_R2.fastq.gz $R1 $R2 --discard-untrimmed; 
done

conda deactivate
exit


###For the samples from LGC company
cd (file position)
for f in *.bz2; do bzcat "$f" | gzip -c >"${f%.*}.gz"; done

