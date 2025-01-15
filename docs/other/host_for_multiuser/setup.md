(set-optinist-config)=
Multi-user (Host setup)
=======================

```{contents}
:depth: 4
```

Follow the steps below to setup `multiuser` mode.

### Clone the Repository
1. In your hosting server, clone the OptiNiSt repository.
    ```bash
    git clone git@github.com:oist/optinist.git -b main
    ```
2. copy config files
    ```bash
    cp -i studio/config/.env.example studio/config/.env
    cp -i studio/config/auth/firebase_config.example.json studio/config/auth/firebase_config.json
    ```

###  Setup Firebase Authentication

#### Create your Firebase Project
1. Go to [https://console.firebase.google.com/](https://console.firebase.google.com/)
2. Click "Add project".
3. Enter your project name, and click "Continue".
4. Google Analytics is optional. You can choose "Enable Google Analytics for this project" or not.
5. After your project is ready, click "Continue".

#### Setup Firebase Authentication
1. Select "Build > Authentication" from the left menu.
2. Select "Get started".
3. Select "Sign-in method" tab.
4. Select "Add new provider" in "Sign-in providers" section.
5. Click "Email/Password" and enable it.
6. Click "Save".

#### Create Admin User for the Project
1. Select "Authentication" from the left menu.
2. Select "Users" tab.
3. Click "Add user" button.
4. Fill the form and click "Add user".
    - Email address: your email address
    - Password: your password

- Created user's "User UID" is required later.

#### Get Firebase Tokens
1. Click setting icon(besides Project Overview), then select "Project settings" from the left menu.
2. Select "General" tab.
3. Select "web app" in "Your apps" section.
4. Enter your app name, and click "Register app".
    - "Also set up Firebase Hosting for this app" is not required.
5. Click continue to console.
6. Set the following values to `studio/config/auth/firebase_config.json`.
    - apiKey
    - authDomain
    - projectId
    - storageBucket
    - messagingSenderId
    - appId
    - (keep databaseURL blank)
7. Select "Service accounts" tab.
8. Click "Generate new private key" in "Firebase Admin SDK" section.
9. Save the downloaded file to `studio/config/auth/firebase_private.json`.

### Setup Database
- Set up your own mysql (or mariadb) server or use docker compose mysql
- Below are the instructions for using mysql with docker compose.

1. Edit configs.
    - Edit studio/config/.env
      - Set `MYSQL_SERVER` to db server host or ip
        - Format: `{DB_HOST}:{DB_PORT}`
        - \*For docker platform, the fixed value `db:3306` is fine.
      - Set `MYSQL_ROOT_PASSWORD` to database root password, which you have decided.
      - Set `MYSQL_DATABASE` to `{YOUR_DATABASE_NAME}`, which you have decided.
      - Set `MYSQL_USER` to `{DB_USER_NAME}`, which you have decided.
      - Set `MYSQL_PASSWORD` to `{DB_USER_PASSWORD}`, which you have decided.
2. Install & run mysql server.
    ```bash
    docker compose -f docker-compose.dev.multiuser.yml up db -d
    ```
    - The database and db_user are automatically generated based on the .env settings.
3. Check connection to mysql server.
    - Connecting via docker command
      ```bash
      docker exec -it {DB_CONTAINER_NAME} mysql -u {DB_USER_NAME} -p {YOUR_DATABASE_NAME}
      mysql> exit
      ```
      - Note: `{DB_CONTAINER_NAME}` is the container name or container ID of the database docker container. (Can be confirmed with `docker ps`)
    - Connect via mysql command (requires mysql-client)
      ```bash
      mysql -h {DB_HOST} --port={DB_PORT} -u {DB_USER_NAME} -p {YOUR_DATABASE_NAME}
      mysql> exit
      ```
    - If a connection to the database server is available, the setup was successful.

### Setup & Run OptiNiSt

#### For Docker Platform

To use multiuser mode with Docker, perform the following steps.

##### Setup Backend

(set-optinist-config)=
###### Set OptiNiSt Config
- Edit `studio/config/.env`
    - Change `SECRET_KEY` to any random string.
    - Change `USE_FIREBASE_TOKEN` to `True`.
    - Change `IS_STANDALONE` to `False`

###### Start Backend (Database is set up on startup)
```bash
docker compose -f docker-compose.dev.multiuser.yml up studio-dev-be -d
```
(insert_initial_data)=
###### Insert Initial Data
```bash
docker exec -it {DB_CONTAINER_NAME} mysql -u {DB_USER_NAME} -p {YOUR_DATABASE_NAME}
```
Make an initial sql entry
```sql
INSERT INTO organization (name) VALUES ('{YOUR_ORGANIZATION_NAME}');
INSERT INTO roles (id, role) VALUES (1, 'admin'), (20, 'operator');
INSERT INTO users (uid, organization_id, name, email, active) VALUES ('{FIREBASE_USER_UID}', 1, '{YOUR_NAME}', '{YOUR_EMAIL}', true);
INSERT INTO user_roles (user_id, role_id) VALUES (1, 1);
```
  - Note on Variables
    - `{FIREBASE_USER_UID}` ... The user uid you created in the previous step ([Create admin user for the project](#create-admin-user-for-the-project)).
    - `{YOUR_ORGANIZATION_NAME}` ... Display name on system (Any text)
    - `{YOUR_NAME}` ... Display name on system (Any text)
    - `{YOUR_EMAIL}` ... Email address corresponding to `{FIREBASE_USER_UID}`

- About Roles
  - Only 2 roles, "admin" and "operator" are supported for now.
    - "admin"
      - can manage other users
    - "operator"
      - general user
  - More information is [here](usage.md).

##### Run OptiNiSt
```bash
docker compose -f docker-compose.dev.multiuser.yml up -d
```

1. Access to `http://{YOUR_HOST}:8000` from your browser.
2. Confirm that you can SingIn with your Firebase Authentication account.

#### For Non-Docker Platforms

Below are the steps for a case using Non-Docker platforms (Windows, Mac, Linux).

##### Setup Backend
- See [OptiNiSt installation guide](../../installation/index.rst) for installing conda and optinist for installing conda and optinist.
- After creating and activating a conda environment for the project, run following commands

###### Set OptiNiSt Config
- Follow the steps to change the `.env` file to multi-user mode, see [Set OptiNiSt config](#set-optinist-config).

###### Setup Database
```bash
cd {OPTINIST_ROOT_PATH}  # root path of repository cloned
alembic upgrade head
```

###### Insert Initial Data
- To setup the database we need to insert some initial data. Follow the procedure in [Insert Initial Data](#insert_initial_data), but remove the initial command `docker exec -it {DB_CONTAINER_NAME}`, as not using Docker. Instead to add initial data and setup database, use:
```bash
mysql -u {DB_USER_NAME} -p {YOUR_DATABASE_NAME}
```

##### Run OptiNiSt
```bash
python main.py
```

1. Access to `http://{YOUR_HOST}:8000` from your browser (most likely `http://localhost:8000`).
2. Confirm that you can SingIn with your Firebase Authentication account.
