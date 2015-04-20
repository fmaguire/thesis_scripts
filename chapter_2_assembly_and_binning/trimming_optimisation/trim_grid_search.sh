#!/bin/bash
#./trim_grid_search.sh PARAM_FILE

usage="./trim_grid_search.sh PARAM_FILE"
if [ "$#" -ne 1 ]
then 
    echo "$usage"
    exit 1
fi

INDEX=/storage/fin/assemblies/trim_illumina_term12/trinity_bulk_reassembly_with_current_trinity/bowtie2_index/bulk_reassembly_trinity
MIN_INS=37
MAX_INS=1161
param_file=$1
summary_file=${param_file}_summary

echo "param_index, sample_id, PE_reads_mapping_concordantly exactly 1 time" > $summary_file

readarray PARAMS < $param_file
#ILLUMINACLIP:/storage/fin/adapters/exeter_sequencing_adaptors.fasta:2 5: 25 30 35 40: 
#LEADING: 5 15 20 25 30 35
#TRAILING: 5 15 20 25 30 35
#SLIDINGWINDOW: 5 10 15 20 : 5 10 15 20 25 30 35

for i in "${!PARAMS[@]}"
do
    for id in Sample*
    do
        R1=${id}/*_R1_sample.fq
        R2=${id}/*_R2_sample.fq
        java -jar /storage/fin/other_tools/Trimmomatic-0.32/trimmomatic-0.32.jar  \
        PE \
        -threads 22 \
        -phred33 \
        $R1 \
        $R2 \
        ${id}/${i}_R1_paired_output.fastq \
        ${id}/${i}_R1_unpaired_output.fastq \
        ${id}/${i}_R2_paired_output.fastq \
        ${id}/${i}_R2_unpaired_output.fastq \
        ${PARAMS[$i]} > /dev/null 2>&1
       
        temp_stats=${id}/${param_file}_${i}.stats
        bowtie2 -q -p 22 -x ${INDEX} \
        -1 ${id}/${i}_R1_paired_output.fastq \
        -2 ${id}/${i}_R2_paired_output.fastq \
        -I ${MIN_INS} -X ${MAX_INS} 2> $temp_stats 1> /dev/null
        
        PE_1_time=$(awk '{print $1}' <(grep "aligned concordantly exactly 1 time" $temp_stats))
        SE_1_time=$(awk '{print $1}' <(grep "aligned discordantly 1 time" $temp_stats))
        Reads_mapped=$(awk '{print $1}' <(grep "reads; of these" $temp_stats))
        

        rm ${id}/${i}_R1_paired_output.fastq \
            ${id}/${i}_R1_unpaired_output.fastq \
            ${id}/${i}_R2_paired_output.fastq \
            ${id}/${i}_R2_unpaired_output.fastq
        
        echo "$i, $id, $Reads_mapped, $PE_1_time, $SE_1_time" 
        echo "$i, $id, $Reads_mapped, $PE_1_time, $SE_1_time" >> $summary_file
    done
done
