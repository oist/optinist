from sqlalchemy.orm import sessionmaker
from sqlmodel import create_engine

from studio.app.common.core.mode import MODE
from studio.app.common.db.config import DATABASE_CONFIG

if not MODE.IS_STANDALONE:
    engine = create_engine(DATABASE_CONFIG.DATABASE_URL, pool_recycle=360, pool_size=20)
    SessionLocal = sessionmaker(
        autocommit=False, autoflush=False, expire_on_commit=False, bind=engine
    )


def get_db():
    if not MODE.IS_STANDALONE:
        try:
            db = SessionLocal()
            yield db
        finally:
            db.close()
