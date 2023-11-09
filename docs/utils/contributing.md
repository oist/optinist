(Contributing)=
# Contributing to OptiNiSt
OptiNiSt welcomes your contributions.
See the following guidelines before submitting Pull Requests.

## Coding Style
### Pre Commit
- run following command before your first commit
  ```
  pre-commit install
  ```
  - Once installed, it automatically checks your code style on every commits.

### Python
- Format all files using [black](https://black.readthedocs.io/en/stable/#).
- Check the code for problems using [flake8](https://pypi.org/project/flake8/).
  - Some excluded rules are contained in `.flake8`.
  - These exclusion are to avoid conflicts between black and flake8.
- Sort your import using [isort](https://github.com/PyCQA/isort).
- These guides are checked by pre-commits, github workflow on submitting Pull Requests.

#### VSCode Settings
- If you are using VSCode, you can use features (like format/lint on save) by installing and activating the extensions in `.vscode/extensions.json`.
- For these extensions, we provide example settings in `.vscode/settings.example.json`
  - To use this,
    - `cp .vscode/settings.example.json .vscode/settings.json`
    - then edit `.vscode/settings.json` following the comment in the file.
      - Uncomment the line `// "flake8.path": ["conda", "run", "-n", "optinist_dev", "python", "-m", "flake8"],`
      - and set `optinist_dev` to your conda environment name.

#### Other editors
- Though automatically checked with pre-commit, you can check your code follow the guide before commit.
- Run the following command to check linting and format your code.
  ```
  make format
  ```
