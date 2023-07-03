from typing import Any, Dict, Optional, Union

from pydantic import BaseModel


class NWBParams(BaseModel):
    session_description: str = "optinist"
    identifier: str = "optinist"
    experiment_description: Optional[str] = None
    device: Union[Dict, Any]
    optical_channel: Union[Dict, Any]
    imaging_plane: Union[Dict, Any]
    image_series: Union[Dict, Any]
    ophys: Union[Dict, Any]
