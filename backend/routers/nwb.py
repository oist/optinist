
from fastapi import APIRouter
import os
import yaml
from .utils import get_dict2nest

router = APIRouter()

nwb_config = {
  'session_description': 'optinist',
  'identifier': 'optinist',
  'experiment_description': 'None',
  'device': {
      'name': 'Microscope device',
      'description': 'Microscope Information',
      'manufacturer': 'Microscope Manufacture'
  },
  'optical_channel': {
      'name': 'OpticalChannel',
      'description': 'optical channel',
      'emission_lambda': 500.1
  },
  'imaging_plane': {
      'name': 'ImagingPlane',
      'description': 'standard',
      'imaging_rate': 30.1,
      'excitation_lambda': 600.1,
      'indicator': 'GCaMap',
      'location': 'V1',
  },
  'image_series': {
      'name': 'TwoPhotonSeries',
      'starting_time': 0,
      'starting_frame': 0,
  },
  'ophys': {
      'plane_segmentation': {
        'children': {
          'name': 'PlaneSegmentation',
          'description': '',
        }
      }
  }
}


@router.get("/nwb")
async def params():
    config = get_dict2nest(nwb_config)
    print(config)
    return config
