from datetime import datetime
from typing import Optional

from sqlalchemy.dialects.mysql import BIGINT
from sqlalchemy.sql.functions import current_timestamp
from sqlmodel import Column, Field, SQLModel, text


class Base(SQLModel):
    id: int = Field(sa_column=Column(BIGINT(unsigned=True), primary_key=True))


class TimestampMixin(SQLModel):
    created_at: Optional[datetime] = Field(
        sa_column_kwargs={"server_default": current_timestamp()},
    )

    updated_at: Optional[datetime] = Field(
        sa_column_kwargs={
            "server_default": text("CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP"),
        },
    )
