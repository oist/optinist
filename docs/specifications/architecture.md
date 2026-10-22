## Architecture
### Frontend
- language: TypeScript
- environment: Node.js (v20)
- framework: React (v18)
- middleware: Redux, Redux Toolkit, React Router

#### Directory Structure
Following is the summary of directory structure.

- frontend
  - src
    - api: handle API requests.
    - components: React components.
    - const: constants.
    - pages: Components set as `element` in `Route` component of react-router.
    - store: Redux store.
      - slice
        - Some
          - SomeActions.ts: create actions like async thunk from api response.
          - SomeSelectors.ts: select state from store.
          - SomeSlice.ts: define reducers.
          - SomeType.ts: define data type.
          - SomeUtils.ts: define utility functions if needed.
    - style: extra css files.
    - utils: common utility functions.
  - `.env`: Application environment variable definition file
    - Depending on usage, the following .env files are available
      - `.env`: Loads for all environments
      - `.env.development`: Used during development (yarn start)
      - `.env.production`: Used during production builds (yarn build)
    - [Learn more](https://create-react-app.dev/docs/adding-custom-environment-variables/#what-other-env-files-can-be-used)

- api, components, store/slice directories are separated by its domain.

### Backend
- language: Python (v3.8)
- framework: FastAPI
- handle workflow by snakemake.

#### Directory Structure
Following is the summary of directory structure.
- studio
  - alembic: files for database migration. database is used in multi-user mode.
  - app
    - `Snakefile`: File read by snakemake. It defines how to run workflow by language like Python.
    - common: common modules(not for specific domain)
      - core: define core functions.
      - dataclass: define dataclasses that are passed between nodes or used as visualize outputs.
      - db: database general configuration.
      - models: database models.
      - routers: define fastapi routers.
      - schemas: define pydantic schemas for api requests and responses.
      - wrappers: functions for algorithm nodes.
    - optinist: modules specific to calcium imaging data processing domain. It has same structure as common.
  - config: configuration files.
    - `.env`: Application environment variable definition file
    - logging:
      - `logging.yaml`: logging config for standalone mode. In standalone mode, logs are sent to streams.
      - `logging.multiuser.yaml`: logging config for multi-user mode. In multi-user mode, logs are sent to file (in `logs` dir). log file rotation settings are included.
      - you can manage log configurations by editing these files.
    - auth: put Firebase authentication secret files.
    - docker: backend Dockerfiles for each environment.
  - test_data: test data for unit tests.
  - tests: unit tests.
  - `__main__.py`: entry point of backend in production.
  - `__main_unit__.py`: core main unit. called by `studio/__main__.py` and `main.py`.
- `main.py`: entry point of backend in development.
- `run_cluster.py`: entry point for CLI run.
