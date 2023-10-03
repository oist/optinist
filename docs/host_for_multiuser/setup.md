Host Optinist yourself for multiuser use
========================================

```{contents}
:depth: 4
```

## Setup multiuser mode

Follow the steps below to setup `multiuser` mode.

### Clone the repository
1. In your hosting server, clone the OptiNiSt repository.
    ```bash
    git clone git@github.com:oist/optinist.git -b main
    ```
2. copy config files
    ```bash
    cp frontend/.env.example frontend/.env
    cp studio/config/.env.example studio/config.env
    cp studio/config/auth/firebase_config.example.json studio/config/auth/firebase_config.json
    ```

### Create your Firebase Project
1. Go to [https://console.firebase.google.com/](https://console.firebase.google.com/)
2. Click "Add project".
3. Enter your project name, and click "Continue".
4. Google Analytics is optional. You can choose "Enable Google Analytics for this project" or not.
5. After your project is ready, click "Continue".

### Setup Firebase Authentication
1. Select "Build > Authentication" from the left menu.
2. Select "Get started".
3. Select "Sign-in method" tab.
4. Select "Add new provider" in "Sign-in providers" section.
5. Click "Email/Password" and enable it.
6. Click "Save".

### Create admin user for the project
1. Select "Authentication" from the left menu.
2. Select "Users" tab.
3. Click "Add user" button.
4. Fill the form and click "Add user".
    - Email address: your email address
    - Password: your password

- Created user's "User UID" is required later.

### Get Firebase tokens
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
9. Put the downloaded file to `studio/config/auth/firebase_private.json`.

### Setup mysql(or mariadb)
1. Install mysql(or mariadb) server.
2. Connect mysql server.
    ```bash
    mysql -u root -p
    ```
3. Create database for your project.
    ```sql
    CREATE DATABASE YOUR_DATABASE_NAME;
    ```
4. Create user for your project.
    ```sql
    CREATE USER 'DB_USER_NAME'@'localhost' IDENTIFIED BY 'DB_USER_PASSWORD';
    ```
5. Grant all privileges to the user for the database.
    ```sql
    GRANT ALL PRIVILEGES ON YOUR_DATABASE_NAME.* TO 'DB_USER_NAME'@'localhost';
    ```

### Set OptiNiSt config
1. Edit `frontend/.env`
    - Change `REACT_APP_IS_STANDALONE` to `false`
2. Edit `studio/config/.env`
    - Change `SECRET_KEY` to any random string.
    - Change `USE_FIREBASE_TOKEN` to `True`.
    - Change `IS_STANDALONE` to `False`
    - Set `MYSQL_SERVER` to your host
    - Set `MYSQL_DATABASE` to {YOUR_DATABASE_NAME}, which you created in the previous step.
    - Set `MYSQL_USER` to {DB_USER_NAME}, which you created in the previous step.
    - Set `MYSQL_PASSWORD` to {DB_USER_PASSWORD}, which you created in the previous step.
    - `MYSQL_ROOT_PASSWORD` can be left commented.

### Setup Frontend
1. Install node.js(version 18)
    - https://nodejs.org/ja
2. Install yarn
    ```bash
    npm install -g yarn
    ```
3. Install frontend requirements
    ```bash
    cd frontend
    yarn install
    ```
4. Build frontend
    ```bash
    yarn build
    ```

### Setup Backend
- See OptiNiSt installation guide.
- After create and activate conda environment for the project, run following commands

1. Install backend requirements
    ```bash
    cd studio
    pip install .
    ```
2. Setup database
    ```bash
    alembic upgrade head
    ```
3. Insert initial data

    ```bash
    mysql -u DB_USER_NAME -p
    ```
    ```sql
    USE YOUR_DATABASE_NAME;
    INSERT INTO organization (name) VALUES ('YOUR_ORGANIZATION_NAME');
    INSERT INTO roles (id, role) VALUES (1, 'admin'), (20, 'operator');
    INSERT INTO users (uid, organization_id, name, email, active, ) VALUES ('USER_UID', 1, 'YOUR_EMAIL', 'YOUR_PASSWORD', 1);
    INSERT INTO user_roles (user_id, role_id) VALUES (1, 1);
    ```
    - USER_UID is the user uid you created in the previous step ([Create admin user for the project](#create-admin-user-for-the-project)).
    - Only 2 roles, "admin" and "operator" are supported for now.
      - "admin"
        - can manage other users
      - "operator"
        - general user

### Run OptiNiSt
```bash
python main.py
```

- Access to `http://{YOUR_HOST}:8000` from your browser.
