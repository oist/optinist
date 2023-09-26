from typing import List

from fastapi import APIRouter, Depends, Query
from sqlmodel import Session, or_

from studio.app.common.core.auth.auth_dependencies import get_current_user
from studio.app.common.db.database import get_db
from studio.app.common.models import User as UserModel
from studio.app.common.schemas.users import User, UserInfo

router = APIRouter(prefix="/users/search", tags=["users/search"])


@router.get(
    "/share_users",
    response_model=List[UserInfo],
    description="""
- Get a list of users with whom content is shared.
- Note: Maximum of 10 responses. (security considerations)
""",
)
def search_share_users(
    keyword: str = Query(description="partial match (user.name or user.email)"),
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    MAX_RESPONSE_COUNT = 10
    users = (
        db.query(UserModel)
        .filter(
            UserModel.active.is_(True),
            UserModel.organization_id == current_user.organization.id,
        )
        .filter(
            or_(
                UserModel.name.like("%{0}%".format(keyword)),
                UserModel.email.like("%{0}%".format(keyword)),
            )
        )
        .order_by(UserModel.id)
        .limit(MAX_RESPONSE_COUNT)
        .all()
    )
    return users
