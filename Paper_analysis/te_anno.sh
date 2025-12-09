#!/bin/bash
##################
# slurm settings #
##################
#SBATCH --job-name te_anno
#SBATCH --qos=vshort
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8000
#SBATCH --array=1-12
#SBATCH --output=/users/asebe/smontgomery/tmp/logs/%x.%a.out
#SBATCH --error=/users/asebe/smontgomery/tmp/logs/%x.%a.err
##################################
# make bash behave more robustly #
##################################
# set -e
# set -u
# set -o pipefail
##Adapted from Shuangyang and Elin to make gtf files for TE annotations
##Made complicated by Sean
# === begin ENVIRONMENT SETUP ===
##Loop through species
if [[ ${SLURM_ARRAY_TASK_ID} == 1 ]]; then
    species=Acas #Acanthamoeba castellani; 14 marks (H3,input)
    chromsize=42019824
elif [[ ${SLURM_ARRAY_TASK_ID} == 2 ]]; then
    species=Atha #Arabidopsis thaliana; 11 marks (H3,input)
    chromsize=119667750
elif [[ ${SLURM_ARRAY_TASK_ID} == 3 ]]; then
    species=Ddis #Dictyostelium discoideum; 14 marks (H3, input)
    chromsize=34134454
elif [[ ${SLURM_ARRAY_TASK_ID} == 4 ]]; then
    species=Ngru #Naegleria gruberii; 10 marks (H3, input)
    chromsize=40964085
elif [[ ${SLURM_ARRAY_TASK_ID} == 5 ]]; then
    species=Nvec #Nematostella vectensis; 14 marks (H3, input)
    chromsize=269418438
elif [[ ${SLURM_ARRAY_TASK_ID} == 6 ]]; then
    species=Ppat #Physcomitrium patens; 11 marks (H3, input)
    chromsize=471852792
elif [[ ${SLURM_ARRAY_TASK_ID} == 7 ]]; then
    species=Spun #Spizellomyces punctatus; 14 marks (H3,input)
    chromsize=24131112
elif [[ ${SLURM_ARRAY_TASK_ID} == 8 ]]; then
    species=Tthe #Tetrahymena thermophila; 13 marks (H3, input)
    chromsize=103014375
elif [[ ${SLURM_ARRAY_TASK_ID} == 9 ]]; then
    species=Bnat #Bigelowiella natans; 14 marks (H3, input)
    chromsize=94701163
elif [[ ${SLURM_ARRAY_TASK_ID} == 10 ]]; then
    species=Gthe #Guillardia theta; 11 marks (H3, input)
    chromsize=87266873
elif [[ ${SLURM_ARRAY_TASK_ID} == 11 ]]; then
    species=Cfra #Creolimax fragrantissima; 10 marks (H3, input)
    chromsize=44821703
elif [[ ${SLURM_ARRAY_TASK_ID} == 12 ]]; then
    species=Scer #Saccharomyces cerevisiae; 11 marks (H3, input)
    chromsize=12157105
fi
#2. set directory containing sample folders
work_folder=/users/asebe/smontgomery/genomes/${species}
#3. set genome and annotation files
gff3file=/users/asebe/smontgomery/genomes/${species}/EDTA2/*TEanno.gff3
genefile=/users/asebe/smontgomery/genomes/${species}/${species}.gene.all.gff3.gz

# === end ENVIRONMENT SETUP ===

# === begin RUN THE PROGRAM ===

cut -f1,2 /users/asebe/xgraubove/genomes/data/${species}_gDNA.fasta.fai > /users/asebe/smontgomery/genomes/${species}/${species}.chromsizes.txt
awk -v OFS='\t' '{print $1,1,$2}' /users/asebe/smontgomery/genomes/${species}/${species}.chromsizes.txt > /users/asebe/smontgomery/genomes/${species}/${species}.circlize.txt

##File nomenclature:
##${species}.<TE/gene>.<subset>.<filetype>
##species is usually the genome the annotation is based on 
##TE/gene indicates if the annotation is TEs, genes, or both (TE, gene, TEgene)
##Subset is what portion of the annotations of TEs/genes is used (see below)
##all = full annotation
##fragment = non-intact TEs
##intact = intact TEs
##LTR = intact LTR TEs
##nonLTR = intact nonLTR TEs
##PCG = protein coding genes
##Filetype is the type of file (gff3, gtf, bed) <shockedPikachu.jpg>

mkdir -p ${work_folder}
cd ${work_folder}

##TE annoations
##Split gff3 file
cp ${gff3file} ${species}.TE.all.gff3
less -S ${species}.TE.all.gff3 | grep '##\|LTR\|LINE\|SINE' > ${species}.TE.RT.gff3
less -S ${species}.TE.all.gff3 | grep -v 'LTR\|LINE\|SINE' > ${species}.TE.DNA.gff3
less -S ${species}.TE.RT.gff3 | grep -v 'LINE\|SINE' | grep -v 'ID=TE_homo_' > ${species}.TE.RT.LTR.intact.gff3
less -S ${species}.TE.RT.gff3 | grep -v 'LINE\|SINE' | grep '##\|ID=TE_homo_' > ${species}.TE.RT.LTR.fragment.gff3
less -S ${species}.TE.RT.gff3 | grep '##\|LINE' > ${species}.TE.RT.LINE.gff3
less -S ${species}.TE.RT.gff3 | grep '##\|SINE' > ${species}.TE.RT.SINE.gff3
less -S ${species}.TE.DNA.gff3 | grep -v 'ID=TE_homo_' > ${species}.TE.DNA.intact.gff3
less -S ${species}.TE.DNA.gff3 | grep '##\|ID=TE_homo_' > ${species}.TE.DNA.fragment.gff3
##Prep and output gtf files
less -S ${species}.TE.RT.LTR.intact.gff3 \
	| egrep 'Copia_LTR_retrotransposon|Gypsy_LTR_retrotransposon|LTR_retrotransposon|Bel_Pao_LTR_retrotransposon' \
	| sed -r 's/Parent=.*Name/Name/' \
	| sed 's/;Sequence_ontology.*//' \
	| sed 's/;/\t/g' \
	| sed 's/Classification=Unknown/Classification=Unknown\/Unknown/' \
	| sed 's/Classification=//' \
	| sed 's/\//\t/' \
	| awk '{print $1"\t"$2"\tTE_RT_LTR_intact\t"$4"\t"$5"\t.\t"$7"\t.\t""gene_id \""$10"\";","transcript_id \""$9"\";","class_id \"RT\";","order_id \""$11"\";","family_id \""$12"\";"}' \
	| sed 's/Name=//' \
	| sed 's/ID=//' \
	> ${species}.TE.RT.LTR.intact.gtf
less -S ${species}.TE.RT.LTR.fragment.gff3 \
	| grep -v '#' \
	| sed 's/;Sequence_ontology.*//' \
	| sed 's/;/\t/g' \
	| sed 's/Classification=Unknown/Classification=Unknown\/Unknown/' \
	| sed 's/Classification=//' \
	| sed 's/\//\t/' \
	| awk '{print $1"\t"$2"\tTE_RT_LTR_fragment\t"$4"\t"$5"\t.\t"$7"\t.\t""gene_id \""$10"\";","transcript_id \""$9"\";","class_id \"RT\";","order_id \""$11"\";","family_id \""$12"\";"}' \
	| sed 's/Name=//' \
	| sed 's/ID=//' \
	> ${species}.TE.RT.LTR.fragment.gtf
less -S ${species}.TE.RT.LINE.gff3 \
	| grep -v '#' \
	| sed 's/;Sequence_ontology.*//' \
	| sed 's/;/\t/g' \
	| sed 's/Classification=Unknown/Classification=Unknown\/Unknown/' \
	| sed 's/Classification=//' \
	| sed 's/\//\t/' \
	| awk -v OFS='\t' '{if ($12=="Penelope") {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,"PLE",$12} else {print $0}}' \
	| awk '{print $1"\t"$2"\tTE_RT_LINE\t"$4"\t"$5"\t.\t"$7"\t.\t""gene_id \""$10"\";","transcript_id \""$9"\";","class_id \"RT\";","order_id \""$11"\";","family_id \""$12"\";"}' \
	| sed 's/Name=//' \
	| sed 's/ID=//' \
	> ${species}.TE.RT.LINE.gtf
less -S ${species}.TE.RT.SINE.gff3 \
	| grep -v '#' \
	| sed 's/;Sequence_ontology.*//' \
	| sed 's/;/\t/g' \
	| sed 's/Classification=Unknown/Classification=Unknown\/Unknown/' \
	| sed 's/Classification=SINE?/Classification=SINE\/Uncertain/' \
	| sed 's/Classification=//' \
	| sed 's/\//\t/' \
	| awk '{print $1"\t"$2"\tTE_RT_SINE\t"$4"\t"$5"\t.\t"$7"\t.\t""gene_id \""$10"\";","transcript_id \""$9"\";","class_id \"RT\";","order_id \""$11"\";","family_id \""$12"\";"}' \
	| sed 's/Name=//' \
	| sed 's/ID=//' \
	> ${species}.TE.RT.SINE.gtf
less -S ${species}.TE.DNA.intact.gff3 \
	| grep -v '#' \
	| sed 's/;Sequence_ontology.*//' \
	| sed 's/;/\t/g' \
	| sed 's/Classification=Unknown/Classification=Unknown\/Unknown/' \
	| sed 's/Classification=//' \
	| sed 's/\//\t/' \
	| awk -v OFS='\t' '{if ($11=="DNA") {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$12,$12} else {print $0}}' \
	| awk '{print $1"\t"$2"\tTE_DNA_intact\t"$4"\t"$5"\t.\t"$7"\t.\t""gene_id \""$10"\";","transcript_id \""$9"\";","class_id \"DNA\";","order_id \""$11"\";","family_id \""$12"\";"}' \
	| sed 's/Name=//' \
	| sed 's/ID=//' \
	| sed 's/order_id \"DTA/order_id \"TIR/' \
	| sed 's/order_id \"DTC/order_id \"TIR/' \
	| sed 's/order_id \"DTH/order_id \"TIR/' \
	| sed 's/order_id \"DTM/order_id \"TIR/' \
	| sed 's/order_id \"DTT/order_id \"TIR/' \
	| sed 's/\"DTA/\"hAT/' \
	| sed 's/\"DTC/\"CACTA/' \
	| sed 's/\"DTH/\"PIF_Harbinger/' \
	| sed 's/\"DTM/\"Mutator/' \
	| sed 's/\"DTT/\"Tc1_Mariner/' \
	> ${species}.TE.DNA.intact.gtf
less -S ${species}.TE.DNA.fragment.gff3 \
	| grep -v '#' \
	| sed 's/;Sequence_ontology.*//' \
	| sed 's/;/\t/g' \
	| sed 's/Classification=Unknown/Classification=Unknown\/Unknown/' \
	| sed 's/Classification=//' \
	| sed 's/\//\t/' \
	| awk -v OFS='\t' '{if ($11=="DNA") {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$12,$12} else {print $0}}' \
	| awk '{print $1"\t"$2"\tTE_DNA_fragment\t"$4"\t"$5"\t.\t"$7"\t.\t""gene_id \""$10"\";","transcript_id \""$9"\";","class_id \"DNA\";","order_id \""$11"\";","family_id \""$12"\";"}' \
	| sed 's/Name=//' \
	| sed 's/ID=//' \
	| sed 's/order_id \"DTA/order_id \"TIR/' \
	| sed 's/order_id \"DTC/order_id \"TIR/' \
	| sed 's/order_id \"DTH/order_id \"TIR/' \
	| sed 's/order_id \"DTM/order_id \"TIR/' \
	| sed 's/order_id \"DTT/order_id \"TIR/' \
	| sed 's/\"DTA/\"hAT/' \
	| sed 's/\"DTC/\"CACTA/' \
	| sed 's/\"DTH/\"PIF_Harbinger/' \
	| sed 's/\"DTM/\"Mutator/' \
	| sed 's/\"DTT/\"Tc1_Mariner/' \
	> ${species}.TE.DNA.fragment.gtf
cat ${species}.TE.RT.LTR.intact.gtf ${species}.TE.RT.LTR.fragment.gtf ${species}.TE.RT.LINE.gtf ${species}.TE.RT.SINE.gtf | sort -k1,1 -k4,4n > ${species}.TE.RT.gtf
cat ${species}.TE.DNA.intact.gtf ${species}.TE.DNA.fragment.gtf | sort -k1,1 -k4,4n > ${species}.TE.DNA.gtf
cat ${species}.TE.RT.LTR.intact.gtf ${species}.TE.DNA.intact.gtf | sort -k1,1 -k4,4n > ${species}.TE.intact.gtf
cat ${species}.TE.RT.gtf ${species}.TE.DNA.gtf | sort -k1,1 -k4,4n > ${species}.TE.all.gtf
##Make bed files from gtf files with length of TE
for f in ${species}.TE.*.gtf; do \
	NAME=`basename $f .gtf`; \
	sed 's/\"//g' $f \
	| awk -v OFS=' ' '{print $1,$4-1,$5,$12,$6,$7,$2,$3,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18"\t"$5-$4}' \
	| sed 's/;//' \
	| sed 's/ /\t/' \
	| sed 's/ /\t/' \
	| sed 's/ /\t/' \
	| sed 's/ /\t/' \
	| sed 's/ /\t/' \
	| sed 's/ /\t/' \
	| sed 's/ /\t/' \
	| sed 's/ /\t/' \
	| sed 's/ /\t/' \
	> ${NAME}.bed;
done


##Gene annotations
##Get full full annotations (w/ tRNA, miRNA, rRNA)
cp ${genefile} ${species}.gene.raw.gtf.gz
zcat ${species}.gene.raw.gtf.gz > ${species}.gene.raw.gtf
zcat ${species}.gene.raw.gtf.gz \
	| grep -v '#' --text \
	| grep -v 'transposable_element' \
	| grep -v 'tRNA_pseudogene' \
	| grep '	gene	' \
	| sed 's/gene_name.*gene_biotype/gene_biotype/g' \
	| sed 's/gene_source.*gene_biotype/gene_biotype/g' \
	| sed 's/\"//g' \
	| sed 's/;//g' \
	| sed 's/lnc_RNA/lncRNA/g' \
	| awk -v OFS='\t' -v sp="${species}" '{if ($9=="gene_id") {print $1,$2,$12,$4,$5,$6,$7,$8,$9" ""\""sp"_"$10"\";"} else {print $1,$2,"protein_coding",$4,$5,$6,$7,$8,"gene_id \""sp"_"$9"\";"}}' \
	| grep -v 'ena	.*RNA	' \
	> ${species}.gene.all.gtf
##Extract exons
zcat ${species}.gene.raw.gtf.gz \
	| grep -v '#' --text \
	| grep -v 'transposable_element' \
	| grep -v 'tRNA_pseudogene' \
	| grep '	exon	' \
	| sed 's/gene_name.*gene_biotype/gene_biotype/g' \
	| sed 's/gene_source.*gene_biotype/gene_biotype/g' \
	| sed 's/\"//g' \
	| sed 's/;//g' \
	| sed 's/lnc_RNA/lncRNA/g' \
	| awk -v OFS='\t' -v sp="${species}" '{if ($9=="gene_id") {print $1,$2,$12,$4,$5,$6,$7,$8,$9" ""\""sp"_"$10"\";"} else {print $1,$2,"protein_coding",$4,$5,$6,$7,$8,"gene_id \""sp"_"$9"\";"}}' \
	| grep -v 'ena	.*RNA	' \
	> ${species}.gene.exon.gtf
##Split gtf
grep -i 'protein_coding' ${species}.gene.all.gtf | awk -v OFS='\t' '{print $1,$2,"PCG",$4,$5,$6,$7,$8,$9" "$10}' > ${species}.gene.PCG.gtf
##Gene gtf to bed
sed 's/\"//g' ${species}.gene.all.gtf | awk -v OFS='\t' '{print $1,$4-1,$5,$10,$6,$7,$2,$3,$8,$9" "$10,$5-$4}' | sed 's/;//' > ${species}.gene.all.bed
##Split bed
grep -i 'protein_coding' ${species}.gene.all.bed | awk -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7,"PCG",$9,$10" "$11" "$12" "$15,$16}' > ${species}.gene.PCG.bed


##Remove temporary files
rm tmp
rm tmp2
rm sed*
rm ${species}.gene.raw.gtf.gz
find ${work_folder} -maxdepth 1 -type f -empty -print -delete
find ${work_folder} -maxdepth 1 -type f -name "*.gff3" -size -350c -delete