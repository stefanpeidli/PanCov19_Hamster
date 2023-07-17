cd data/

# LiaoZhang120520
mkdir -p LiaoZhang120520
cd LiaoZhang120520
wget -O GSE145926_RAW.tar 'http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE145926&format=file'
tar -xf GSE145926_RAW.tar
rm GSE145926_RAW.tar
cd ..

# MelmsIzar290421
# https://singlecell.broadinstitute.org/single_cell/study/SCP1219/
mkdir -p MelmsIzar290421
cd MelmsIzar290421
# curl "https://singlecell.broadinstitute.org/single_cell/api/v1/bulk_download/generate_curl_config?accessions=SCP1219&auth_code=3KMlHbTr&directory=all"  -o cfg.txt; curl -K cfg.txt && rm cfg.txt
wget -O GSE171524_RAW.tar 'http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE171524&format=file'
tar -xf GSE171524_RAW.tar
rm GSE171524_RAW.tar
curl -O https://ftp.ncbi.nlm.nih.gov/geo/series/GSE171nnn/GSE171524/suppl/GSE171524_lung_metaData.txt.gz
curl -O https://ftp.ncbi.nlm.nih.gov/geo/series/GSE171nnn/GSE171524/suppl/GSE171524_processed_data.csv.gz
cd ..

# Get ChuaEils from here:
# https://figshare.com/articles/dataset/COVID-19_severity_correlates_with_airway_epithelium-immune_cell_interactions_identified_by_single-cell_analysis/12436517