# optinist <img src="docs/_static/optinist.png" width="250" title="optinist" alt="optinist" align="right" vspace = "50">

<p align="center">
    <a>
       <img src="https://img.shields.io/badge/-Python-F9DC3E.svg?logo=python&style=flat">
    </a>
    <a>
      <img src="https://img.shields.io/badge/-TypeScript-007ACC.svg?logo=typescript&style=flat&logoColor=white">
    </a>
    <a href="https://pypi.org/project/optinist/">
        <img alt="PYPI" src="https://static.pepy.tech/personalized-badge/optinist?period=total&units=international_system&left_color=black&right_color=blue&left_text=Downloads(PYPI)">
    </a>
    <a href="https://pypi.org/project/optinist/">
        <img alt="PYPI" src="https://static.pepy.tech/personalized-badge/optinist?period=week&units=international_system&left_color=black&right_color=blue&left_text=Downloads/week(PYPI)">
    </a>
    <a href="https://github.com/oist/optinist">
      <img alt="" src="https://badge.fury.io/py/optinist.svg">
    </a>
    <a href="https://github.com/oist/optinist">
      <img alt="" src="https://img.shields.io/github/repo-size/oist/optinist">
    </a>
    <a href="https://github.com/oist/optinist">
      <img alt="" src="https://img.shields.io/github/stars/oist/optinist?style=social">
    </a>
    <a href="https://github.com/oist/optinist">
      <img alt="" src="https://img.shields.io/github/forks/oist/optinist?style=social">
    </a>
</p>

OptiNiSt(Optical Neuroimage Studio) helps researchers try multiple data analysis methods, visualize the results, and construct the data analysis pipelines easily and quickly. OptiNiSt's data-saving format follows NWB standards.

OptiNiSt also supports reproducibility of scientific research, standardization of analysis protocols, and developments of novel analysis tools as plug-in.

## Key Features
### :beginner: Easy-To-Create Workflow
- **zero-knowledge of coding**: OptiNiSt allows you to create analysis pipelines easily on the GUI.

### :zap: Visualizing analysis results
- **quick visualization**: OptiNiSt supports you visualize the analysis results by plotly.

### :rocket: Managing Workflows
- **recording and reproducing**: OptiNiSt records and reproduces the workflow pipelines easily.


## Installation
Need anaconda or miniconda environment.
```
conda create -n optinist python=3.8
conda activate optinist
```

Install from pip.
```
pip install optinist
```

launch.
```
run_optinist
```

## Documentation
https://optinist.readthedocs.io/en/latest/


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


<!-- ## Citing the Project
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
``` -->
