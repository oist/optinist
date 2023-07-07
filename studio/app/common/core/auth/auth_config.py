from pydantic import BaseSettings, Field

from studio.app.dir_path import DIRPATH


class AuthConfig(BaseSettings):
    REFRESH_TOKEN_EXPIRE_MINUTES: int = Field(
        default=60 * 24 * 1, env="REFRESH_TOKEN_EXPIRE_MINUTES"
    )
    SECRET_KEY: str = Field(default="123456", env="SECRET_KEY")
    USE_FIREBASE_TOKEN: bool = Field(default=True, env="USE_FIREBASE_TOKEN")
    ALGORITHM = "HS256"

    class Config:
        env_file = f"{DIRPATH.CONFIG_DIR}/.env"
        env_file_encoding = "utf-8"


AUTH_CONFIG = AuthConfig()
