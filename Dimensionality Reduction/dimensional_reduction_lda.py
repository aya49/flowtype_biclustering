import time

import numpy as np
from sklearn.decomposition import LatentDirichletAllocation
from util import read_csv, save_csv_data
import csv

target = "flowcap"
if target == 'flowcap':
    data_dir = '/Local-Scratch/CMPT-884-workspace/flowCAP-II/result/feat/file-cell-countAdj.csv'
    # data, sample_number_list = read_csv(data_dir)
elif target == 'IMPC':
    data_dir = '/Local-Scratch/CMPT-884-workspace/IMPC/result/P1/Sanger_SPLEEN/feat/file-cell-countAdj.csv'
    # data, sample_number_list = read_csv(data_dir)
elif target == 'IMPCPEER':
    data_dir = '/Local-Scratch/CMPT-884-workspace/IMPC/result/P1/Sanger_SPLEEN/feat/file-cell-countAdjPEER.csv'
else:
    raise ValueError("wrong target")
print (target)

def apply_lda():
    start_time = time.time()
    data_array_all, sample_number_list = read_csv(data_dir)
    shape = data_array_all.shape
    feature_num = shape[1]
    data_num = shape[0]

    lda = LatentDirichletAllocation(n_components=feature_num / 8)
    lda.fit(data_array_all)
    dimension_reduction = lda.transform(data_array_all)

    save_dir = './dr_results/{1}-lda-{0}.csv'.format(feature_num / 8, target)
    save_csv_data(dir=save_dir, data=dimension_reduction,
                  sample_number_list=sample_number_list,
                  name='lda')
    # return dimension_reduction
    print("---lda for {0} {1} seconds ---".format(target, time.time() - start_time))


if __name__ == '__main__':
    apply_lda()
    print 'still working'
