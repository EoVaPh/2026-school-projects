import os

os.system('mmseqs easy-cluster seqres_all_new.txt cluster_results cluster_tmp --min-seq-id 0.5 -c 0.5 --cov-mode 0 --threads 8')
