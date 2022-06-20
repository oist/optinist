## Implemented Analysis
### Image analysis
  - [CaImAn](https://github.com/flatironinstitute/CaImAn)
  - [suite2p](https://github.com/MouseLand/suite2p)

### PostProcess

#### basic neural analysis
- Event Trigger Average

#### dimension reduction
- [Canonical Correlation Analysis](https://scikit-learn.org/stable/modules/generated/sklearn.cross_decomposition.CCA.html)
- [Principal Component Analysis](https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html?highlight=pca#sklearn.decomposition.PCA)
- [T-distributed Stochastic Neighbor Embedding](https://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html?highlight=tsne#sklearn.manifold.TSNE)

#### neural population analysis
- Correlation
- [Cross Correlation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.correlate.html)
- [Granger causality test](https://www.statsmodels.org/dev/generated/statsmodels.tsa.stattools.grangercausalitytests.html)

#### neural decoding
- [Generalized Linear Models](https://www.statsmodels.org/stable/glm.html)
- [Linear Discriminant Analysis](https://scikit-learn.org/stable/modules/generated/sklearn.discriminant_analysis.LinearDiscriminantAnalysis.html)
- [Support Vector Machines](https://scikit-learn.org/stable/modules/svm.html#svm)