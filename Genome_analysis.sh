path=$1
files=$(ls $path)
for filename in $files
do 
	[ -d $filename ] && echo ${filename} >> samplename.txt
done
while read line
do
i=$line

#lumpy_sv
cd "$i"
source activate py27
bwa mem -R "@RG\tID:11\tSM:11-10\tLB:lib" /home/dqlab/S288C-2u.fa "$i"_1.fq.gz "$i"_2.fq.gz > "$i".sam
samblaster --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20 < "$i".sam | samtools view -S -b - > "$i".bam
samtools view -b -F 1294 "$i".bam > "$i".discordants.unsorted.bam
samtools view -h "$i".bam|/home/dqlab/anaconda3/envs/py27/bin/extractSplitReads_BwaMem -i stdin|samtools view -Sb - > "$i".splitters.unsorted.bam
samtools sort "$i".discordants.unsorted.bam > "$i".discordants.sorted
samtools sort "$i".splitters.unsorted.bam > "$i".splitters.sorted
samtools sort "$i".bam > "$i".sorted
lumpyexpress -B "$i".sorted -S "$i".splitters.sorted -D "$i".discordants.sorted -o "$i".vcf
rm "$i".sam "$i".bam "$i".discordants.sorted "$i".splitters.sorted "$i".discordants.unsorted.bam "$i".splitters.unsorted.bam
conda deactivate

#Filter lumpy data
bedtools intersect -a "$i".vcf -b /home/dqlab/WYlumpy.bed -v > filtered1_"$i".vcf
cat /home/dqlab/header_lines.vcf filtered1_"$i".vcf > filtered2_"$i".vcf
bcftools view -i 'INFO/SR>=8' filtered2_"$i".vcf | grep -v '^chr=mt' '^chr=2micro' > temp_filtered_"$i".vcf 
awk -v filename="$i" '{
    if ($0 !~ /^#/) print $0 "\t" filename; 
    else print $0
  }' temp_filtered_"$i".vcf > temp_with_headers_"$i".vcf

grep -v '^#' temp_with_headers_"$i".vcf > filtered_"$i".vcf
rm filtered1_"$i".vcf filtered2_"$i".vcf temp_filtered_"$i".vcf temp_with_headers_"$i".vcf

#SNVs and InDels calling
samtools mpileup -f /home/dqlab/S288C-2u.fa "$i".sorted>"$i".mp
java -jar /home/dqlab/VarScan.v2.3.9.jar mpileup2snp "$i".mp --min-avg-qual 20 --min-reads2 15 --min-coverage 20 --min-var-freq 0.25 --output-vcf 1 >"$i"SnptoS288c
java -jar /home/dqlab/VarScan.v2.3.9.jar mpileup2indel "$i".mp --min-avg-qual 20 --min-reads2 15 --min-coverage 20 --min-var-freq 0.25 --output-vcf 1 >"$i"InDeltoS288c

#Filter SNVs data
java -jar /home/dqlab/VarScan.v2.3.9.jar compare /home/dqlab/SNPforFilter-FL-4.vcf "$i"SnptoS288c unique2 "$i"snpcompare
awk '{print $1,$2,$2+1,$4,$5,$10}' "$i"snpcompare >"$i"snpc
cat "$i"snpc|tr " " "\t">"$i"cc
bedtools intersect -v -a "$i"cc -b /home/dqlab/RepeatRegions.txt>"$i"snpfilter
rm "$i"cc "$i"snpc
cat "$i"snpfilter |tr ":" "\t" >"$i"snpcf
rm "$i"snpfilter
cat "$i"snpcf|tr "%" " ">"$i"snpcff
awk '$12>25 {print $1,$2,$3,$4,$5,$12,$13,$14,$15,$16,$17,$18,$19}' "$i"snpcff > "$i"snpafterfilter
rm "$i"snpcf "$i"snpcff *compare
#Filter InDels data
java -jar /home/dqlab/VarScan.v2.3.9.jar compare /home/dqlab/InDelsforFilter-FL-4.vcf "$i"InDeltoS288c unique2 "$i"indelcompare
awk '{print $1,$2,$2+1,$4,$5,$10}' "$i"indelcompare >"$i"indelc
cat "$i"indelc|tr " " "\t">"$i"indelcc
bedtools intersect -v -a "$i"indelcc -b /home/dqlab/RepeatRegions.txt>"$i"indelfilter
rm "$i"indelcc "$i"indelc
cat "$i"indelfilter |tr ":" "\t" >"$i"indelcf
rm "$i"indelfilter
cat "$i"indelcf|tr "%" " ">"$i"indelcff
awk '$12>25''{print $1,$2,$3,$4,$5,$12,$13,$14,$15,$16,$17,$18,$19}' "$i"indelcff>"$i"indelafterfilter
rm "$i"indelcff "$i"indelcf *compare

awk -v a=$i '{print $1,$2,$4,$5,$6,$7,$8,$10,$11,$12,$13,a}' "$i"snpafterfilter>"$i"snpmark
awk -v a=$i '{print $1,$2,$4,$5,$6,$7,$8,$10,$11,$12,$13,a}' "$i"indelafterfilter>"$i"indelmark

#SNP data calling (W303-1A and YJM789)
java -jar /home/dqlab/VarScan.v2.3.9.jar compare //home/dqlab/SNPforLOH2023.txt "$i".mp intersect "$i".snp
cut -f 5 "$i".snp>SNPsingle
# SNPforLOH2023 is a file listing all SNPs between YJM789 and W303-1A
java -jar /home/dqlab/VarScan.v2.3.9.jar compare "$i".mp /home/dqlab/SNPforLOH2023.txt intersect "$i".mpintersect
cut -f 1,2,3,4,5 "$i".mpintersect>"$i".mpcut
#"$fold"forLOH will be used to generate LOH pictures using R program
paste "$i".mpcut SNPsingle>"$i"forLOH
rm "$i".snp "$i".mp SNPsingle "$i".mpcut "$i".mpintersect *toS288c *afterfilter

#LOH mapping
cp /home/dqlab/Scripts/LOHpicture.R .
Rscript LOHpicture.R
#Recombination annotation
python3 /home/dqlab/Scripts/auto_annotate2.py --file_dir . --output_dir .
#Mitochondrial and 2micor plasmid sequencing depth
cp /home/dqlab/Scripts/calculate_coverage.py .
python3 calculate_coverage.py
#W303-1A and YJM789 sequencing depth
cp /home/dqlab/Scripts/rDNA-CUP1-CNV.R .
Rscript rDNA-CUP1-CNV.R
#Average 500bp, step size 200bp sequencing depth
cp /home/dqlab/Scripts/readdepth500.py .
python3 readdepth500.py . 500 200 .
#Average 500bp, step size 200bp sequencing depth mapping
cp /home/dqlab/Scripts/readdepth2.R .
Rscript readdepth2.R
#Average 100bp, step size 100bp sequencing depth
cp /home/dqlab/Scripts/readdepth100.py .
python3 readdepth100.py . 100 100 .
cp /home/dqlab/Scripts/readdepth_filtered.py .
#Leave coordinates with sequencing depth 0
python3 readdepth_filtered.py
rm LOHpicture.R readdepth500.py readdepth2.R readdepth100.py readdepth_filtered.py rDNA-CUP1-CNV.R calculate_coverage.py
cd ..
done < samplename.txt

