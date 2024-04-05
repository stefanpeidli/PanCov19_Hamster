# Code for integration of Gold and Dwarf timecourse SCoV2 hamsters and comparison with public human scRNA-seq data
### Before you run the code yourself:
- Make sure you have conda installed (I suggest using [Miniconda](https://docs.conda.io/en/latest/miniconda.html)).
- Go to the config/config.yaml file and change the paths to the data files to your own paths.
- Make sure you have the required packages installed (see config/ environment yaml files).

### To run the code:
- First, go into the snakemaker folder and run the following command to run the snakemake pipeline:
`snakemake --jobs 10 -k --use-conda `. Make sure to allocate a gpu since we use scVI for integration.

### Associated [Preprint on bioRxiv](https://doi.org/10.1101/2023.08.25.551434)
Peidli, S., Nouailles, G., Wyler, E., Adler, J. M., Kunder, S., Voss, A., ... & Praktiknjo, S. D. (2023). Single-cell-resolved interspecies comparison identifies a shared inflammatory axis and a dominant neutrophil-endothelial program in severe COVID-19. bioRxiv, 2023-08.
