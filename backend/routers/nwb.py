
from fastapi import APIRouter
import os
import yaml

router = APIRouter()

@router.get("/nwb")
async def params():
    config = {
      'session_description': 'optinist',
      'identifier': 'optinist',
      'device': {
        'children': {
          'name': 'Microscope device',
          'description': 'Microscope Information',
          'manufacturer': 'Microscope Manufacture'
        }
      },
      'optical channel': {
        'children': {
          'name': 'OpticalChannel',
          'description': 'optical channel',
          'emission_lambda': 500
        }
      },
      'imaging plane': {
        'children': {
          'name': 'ImagingPlane',
          'description': 'standard',
          'imaging_rate': 30,
          'excitation_lambda': 600,
          'indicator': 'GCaMap',
          'location': 'V1',
        }
      },
      'image series': {
        'children': {
          'name': 'TwoPhotonSeries',
          'starting_time': 0,
          'starting_frame': 0,
        }
      },
      'ophys': {
        'children': {
          'plane_segmentation': {
            'children': {
              'name': 'PlaneSegmentation',
              'description': '',
            }
          }
        }
      }
    }
    
    return config
