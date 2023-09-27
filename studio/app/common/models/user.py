from datetime import datetime
from typing import Dict, List, Optional

from sqlalchemy.dialects.mysql import BIGINT
from sqlalchemy.sql.functions import current_timestamp
from sqlmodel import (
    JSON,
    Column,
    Field,
    ForeignKey,
    Relationship,
    String,
    UniqueConstraint,
)

from studio.app.common.models.base import Base, TimestampMixin
from studio.app.common.models.workspace import WorkspacesShareUser


class UserRole(Base, table=True):
    __tablename__ = "user_roles"
    __table_args__ = (UniqueConstraint("user_id", "role_id", name="idx_user_id"),)

    user_id: int = Field(
        sa_column=Column(BIGINT(unsigned=True), ForeignKey("users.id"), nullable=False),
    )
    role_id: int = Field(
        sa_column=Column(BIGINT(unsigned=True), ForeignKey("roles.id"), nullable=False),
    )
    created_at: Optional[datetime] = Field(
        sa_column_kwargs={"server_default": current_timestamp()},
    )


class User(Base, TimestampMixin, table=True):
    __tablename__ = "users"
    __table_args__ = (UniqueConstraint("uid", name="idx_uid"),)

    organization_id: int = Field(
        sa_column=Column(
            BIGINT(unsigned=True), ForeignKey("organization.id"), nullable=False
        ),
    )
    uid: str = Field(sa_column=Column(String(100), nullable=False))
    name: str = Field(sa_column=Column(String(100), nullable=False))
    email: str = Field(sa_column=Column(String(255), nullable=False))
    attributes: Optional[Dict] = Field(default={}, sa_column=Column(JSON))
    active: bool = Field(nullable=False)

    workspace: List["Workspace"] = Relationship(  # noqa: F821
        back_populates="user", sa_relationship_kwargs={"uselist": True}
    )
    organization: "Organization" = Relationship(
        back_populates="users", sa_relationship_kwargs={"uselist": False}
    )
    role: Optional["Role"] = Relationship(
        back_populates="user",
        link_model=UserRole,
        sa_relationship_kwargs=dict(uselist=False),
    )
    workspace_share: List["Workspace"] = Relationship(  # noqa: F821
        back_populates="user", link_model=WorkspacesShareUser
    )


class Role(Base, table=True):
    __tablename__ = "roles"

    role: str = Field(sa_column=Column(String(100), nullable=False))
    created_at: Optional[datetime] = Field(
        sa_column_kwargs={"server_default": current_timestamp()},
    )
    user: List["User"] = Relationship(back_populates="role", link_model=UserRole)


class Organization(Base, table=True):
    __tablename__ = "organization"

    name: str = Field(sa_column=Column(String(100), nullable=False))
    created_at: Optional[datetime] = Field(
        sa_column_kwargs={"server_default": current_timestamp()},
    )

    users: List["User"] = Relationship(
        back_populates="organization", sa_relationship_kwargs={"uselist": True}
    )
