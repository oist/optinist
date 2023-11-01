from typing import List

from pydantic import BaseModel


class DeleteItem(BaseModel):
    uidList: List


class RenameItem(BaseModel):
    new_name: str
