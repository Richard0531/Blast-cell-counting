module load R/4.3.2
module load cellpose/3.0.7
module load cellprofiler/4.2.5

SAMPLE=$1
INPUT=$2
OUTPUT=$3

cp /home/yhsiao/pharmaco_imaging/Cellprofiler/blast_output/nn_0515.model $OUTPUT
cd $INPUT
#bash /home/yhsiao/pharmaco_imaging/Cellprofiler/file_rename.sh $SAMPLE
python /home/yhsiao/pharmaco_imaging/Cellprofiler/cellpose_seg_pipeline.py $INPUT
cellprofiler -c -r -p /home/yhsiao/pharmaco_imaging/Cellprofiler/blast/blast_seg_classified_0519_pipeline_addImages.cppipe -i $INPUT -o $OUTPUT
#Rscript /home/yhsiao/pharmaco_imaging/Cellprofiler/blast_counting_umap.R $SAMPLE $OUTPUT










































