## Pull Request
- GitHub Pull Request page
  - [https://github.com/oist/optinist/pulls](https://github.com/oist/optinist/pulls)

### Pre Commit
- run following command before your first commit
  ```
  pre-commit install
  ```
  - Once installed, it automatically checks your coding style on every commits.

### Branch rules
- You can submit Pull Request by pushing new branch.
  - Make sure the base branch is `develop-main`, and PR is to `develop-main`.
  - You can't push to the `develop-main` branch directly, the branch is protected.
  - Make sure new branch name is in following format (`xxx` is the name of the feature or bug you are working on.).
    - `feature/xxx`
    - `fix/xxx`
