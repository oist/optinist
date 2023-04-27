# Contributing to OptiNiSt
OptiNiSt welcomes your contributions.
See the following guidelines before submitting Pull Requests.

## Coding Style
### Python
- Format all files using [black](https://black.readthedocs.io/en/stable/#).
- Check the code for problems using [flake8](https://pypi.org/project/flake8/).
  - Some excluded rules are contained in `setup.cfg`.
  - These exclusion are to avoid conflicts between black and flake8.
- Sort your import using [isort](https://github.com/PyCQA/isort).
- These guides are checked by workflows on submitting Pull Requests.
- If you are using VSCode, you can follow the above guidelines by installing and activating the extensions in `.vscode/extensions.json`.
