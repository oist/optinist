from contextlib import contextmanager

from sqlalchemy.orm import sessionmaker
from sqlmodel import create_engine

from studio.app.common.db.config import DATABASE_CONFIG

engine = create_engine(
    DATABASE_CONFIG.DATABASE_URL, pool_recycle=360, pool_size=DATABASE_CONFIG.POOL_SIZE
)

SessionLocal = sessionmaker(
    autocommit=False, autoflush=False, expire_on_commit=False, bind=engine
)


@contextmanager
def session_scope():
    session = SessionLocal()
    try:
        yield session
        session.commit()
    except:  # noqa
        session.rollback()
        raise
    finally:
        session.close()


def get_db():
    try:
        db = SessionLocal()
        yield db
    finally:
        db.close()
