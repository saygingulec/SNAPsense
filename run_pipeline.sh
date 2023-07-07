module load matlab/2022a

for CELL_NAME in data/*; do
  sbatch -t 1-0 --mem 16G --output="slurm-%j-${CELL_NAME##*/}.out" --wrap="matlab -batch \"SNAPsense_pipeline(\'${CELL_NAME}\')\""
done

