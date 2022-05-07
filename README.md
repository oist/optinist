# optinist <img src="docs/_static/optinist.png" width="250" title="optinist" alt="optinist" align="right" vspace = "50">


OptiNiSt(Optical Neuroimage Studio) that helps researchers to try multiple data analysis methods, to visualize the results, and to construct data analysis pipelines. Data formats follow NWB standard.


OpitNiSt will allow not only neural data experts but also non-experts to try advanced analysis methods by combining neuroimage and behavioral data.
It will help functional data analysis and pipeline building first by Brain/MINDS researchers but also by neuroscientists and data scientists world wide.
Availability of not only data but also the tool for analysis on Brian/MINDS Data Portal will facilitate wider usage of marmoset functional brain data.
OptiNiSt will also help reproducibility of research, standardization of analysis protocols, and developments of novel analysis tools as plug-in.

## Key Features
### :beginner: Easy-To-Use Workflow
- **zero-knowledge of coding**: OptiNiSt can make analysis pipelines through connecting nodes and run on GUI.

### :zap: Visualize analysis result
- **visualize analysis result**: OptiNiSt also visualizes neural data analysis result by plotly.

### :rocket: Record Workflow
- **reproduce past workflow**: OptiNiSt records and reproduces past workflow pipeline.

## Documentation
- [Tutorial](docs/utils/tutorial.md)
- [How to use gui workflow](docs/gui/workflow.md)
- [How to use gui visualize](docs/gui/visualize.md)


## Install
How to install and use optinist
- [For Windows](docs/installation/windows.md)
- [For Mac](docs/installation/mac.md)
- [For Linux](docs/installation/linux.md)
- [For Docker](docs/installation/docker.md)


## Implemented Analysis
- Image analysis
  - [CaImAn](https://github.com/flatironinstitute/CaImAn)
  - [suite2p](https://github.com/MouseLand/suite2p)
- PostProcess
  - basic neural analysis
    - Event Trigger Average
    - Cell Grouping
  - dimension reduction
    - [Canonical Correlation Analysis](https://scikit-learn.org/stable/modules/generated/sklearn.cross_decomposition.CCA.html)
    - [Principal Component Analysis](https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html?highlight=pca#sklearn.decomposition.PCA)
    - [T-distributed Stochastic Neighbor Embedding](https://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html?highlight=tsne#sklearn.manifold.TSNE)
  - neural population analysis
    - Correlation
    - [Cross Correlation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.correlate.html)
    - [Granger causality test](https://www.statsmodels.org/dev/generated/statsmodels.tsa.stattools.grangercausalitytests.html)
  - neural decoding
    - [Generalized Linear Models](https://www.statsmodels.org/stable/glm.html)
    - [Linear Discriminant Analysis](https://scikit-learn.org/stable/modules/generated/sklearn.discriminant_analysis.LinearDiscriminantAnalysis.html)
    - [Support Vector Machines](https://scikit-learn.org/stable/modules/svm.html#svm)


## Using GUI
### Workflow
- OptiNiSt can make analysis pipelines through connecting nodes and run on GUI. There are many analysis flow combinations. It selects nodes and connect edges to create pipelne flow.
<p align="center">
  <img width="400px" src="docs/_static/workflow/whole.png" alt="workflow" />
</p>



### Visualize
- OptiNiSt also visualizes neural data analysis results by plotly. It supports image, roi, timeseries, heatmap and so on.
<p align="center">
  <img width="400px" src="docs/_static/visualize/whole.png" alt="visualize" />
</p>

### Record Workflow
- OptiNiSt records and reproduces past workflow pipelines. It can download results as nwb format.
<p align="center">
  <img width="400px" src="docs/_static/record/whole.png" alt="record" />
</p>



## Contributors
### Proposers
- Kenji Doya, [OIST Neural Computation Unit](https://groups.oist.jp/ncu)
- Yukako Yamane, [OIST Neural Computation Unit](https://groups.oist.jp/ncu)

### Main Developers
- [Shogo Akiyama](https://github.com/ShogoAkiyama)
- [Yoshifumi Takeshima](https://github.com/Yoshifumi14)

### Support Developers
- [Tatsuya Tanabe](https://github.com/ttya16)
- [Yosuke Kaneko](https://github.com/toto-maru)
- [Syuya Saeki](https://github.com/hiiaka)


## Citing the Project
To cite this repository in publications:
```
@misc{OptiNiSt,
  author = {name},
  title = {title},
  year = {2022},
  publisher = {},
  journal = {},
  howpublished = {},
}
```
