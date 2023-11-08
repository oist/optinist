## Test
- We have unit tests for both frontend and backend. They are automatically run by github workflow on submitting Pull Requests.
- You can also run them locally.
  - Frontend and Backend (Requires docker)
    ```
    make test_run
    ```
  - Frontend only
    ```
    yarn test
    # or yarn test:ci
    ```
  - Backend only
    ```
    pytest
    ```
