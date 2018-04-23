import numpy as np
from fabia import fabia
from fabia import generate_dataset

FB = fabia.FabiaBiclustering(n_clusters=3, thresZ=0.7)
Xn, X, zc, lc = generate_dataset.make_fabia_biclusters(shape=[6, 7], n_clusters=3, sample_range=[0, 10],
                                                       feature_range=[0, 10],
                                                       bg_noise=0.0001, z_bg_noise=0.0001, z_mean=10, z_sd=0.0001,
                                                       l_bg_noise=0.01, l_mean=0.01, l_sd=0.01)
FB.fit(X)
a = FB.biclusters_  # Returns the ``rows_`` and ``columns_`` members.
b = FB.L_  # Weights of each feature in each bicluster (factor loadings).
c = FB.Z_  # Weights of each sample in each bicluster (factor weights).
# FB.n_clusters
indices = FB.get_indices(0)
print indices
print 'working'
