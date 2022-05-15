# optinist <img src="docs/_static/optinist.png" width="250" title="optinist" alt="optinist" align="right" vspace = "50">


OptiNiSt(Optical Neuroimage Studio) helps researchers try multiple data analysis methods, visualize the results, and construct the data analysis pipelines easily and quickly. OptiNiSt's data-saving format follows NWB standards.

OptiNiSt also supports reproducibility of scientific research, standardization of analysis protocols, and developments of novel analysis tools as plug-in.

## Key Features
### :beginner: Easy-To-Create Workflow
- **zero-knowledge of coding**: OptiNiSt allows you to create analysis pipelines easily on the GUI.

### :zap: Visualizing analysis results
- **quick visualization**: OptiNiSt supports you visualize the analysis results by plotly.

### :rocket: Managing Workflows
- **recording and reproducing**: OptiNiSt records and reproduces the workflow pipelines easily.

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
- OptiNiSt allows you to make your analysis pipelines by graph style using nodes and edges on GUI. Parameters for each analysis are easily changeable. 
<p align="center">
  <img width="400px" src="docs/_static/workflow/whole.png" alt="workflow" />
</p>



### Visualize
- OptiNiSt allows you to visualize the analysis results with one click by plotly. It supports a variety of plotting styles.
<p align="center">
  <img width="400px" src="docs/_static/visualize/whole.png" alt="visualize" />
</p>

### Record
- OptiNiSt supports you in recording and reproducing workflow pipelines in an organized manner. 
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
