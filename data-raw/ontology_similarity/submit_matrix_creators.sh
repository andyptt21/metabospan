for i in {1..11}
do
    echo "Submitted create_RaMP_pathway_overlap_matrix_$i.sh"
    sbatch create_RaMP_pathway_overlap_matrix_$i.sh
done
