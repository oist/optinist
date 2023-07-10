import { RootState, store } from 'store/store'

import { selectRunPostData } from '../slice/Run/RunSelectors'

describe('RunSelectors', () => {
  const initialRootState = store.getState()

  test('selectRunPostData', () => {
    expect(selectRunPostData(targetState)).toEqual(expectRunPostData)
  })

  test('selectRunPostData when forceRunList exists', () => {
    expect(selectRunPostData(targetStateForceRunList)).toEqual(
      expectRunPostDataForceRunList,
    )
  })

  const targetState = {
    ...initialRootState,
    algorithmNode: {
      suite2p_roi_16ffq69kkc: {
        functionPath: 'suite2p/suite2p_roi',
        name: 'suite2p_roi',
        params: {
          tau: { type: 'child', value: 1, path: 'tau' },
          fs: { type: 'child', value: 10, path: 'fs' },
          soma_crop: { type: 'child', value: true, path: 'soma_crop' },
          high_pass: { type: 'child', value: 100, path: 'high_pass' },
          sparse_mode: { type: 'child', value: true, path: 'sparse_mode' },
          max_overlap: { type: 'child', value: 0.75, path: 'max_overlap' },
          nbinned: { type: 'child', value: 5000, path: 'nbinned' },
          spatial_scale: { type: 'child', value: 0, path: 'spatial_scale' },
          threshold_scaling: {
            type: 'child',
            value: 1,
            path: 'threshold_scaling',
          },
          max_iterations: { type: 'child', value: 20, path: 'max_iterations' },
          spatial_hp_detect: {
            type: 'child',
            value: 25,
            path: 'spatial_hp_detect',
          },
          preclassify: { type: 'child', value: 0, path: 'preclassify' },
          allow_overlap: { type: 'child', value: false, path: 'allow_overlap' },
          inner_neuropil_radius: {
            type: 'child',
            value: 2,
            path: 'inner_neuropil_radius',
          },
          min_neuropil_pixels: {
            type: 'child',
            value: 350,
            path: 'min_neuropil_pixels',
          },
          neucoeff: { type: 'child', value: 0.7, path: 'neucoeff' },
        },
        isUpdated: false,
      },
      suite2p_file_convert_m58owcejm0: {
        functionPath: 'suite2p/suite2p_file_convert',
        name: 'suite2p_file_convert',
        params: {
          nplanes: { type: 'child', value: 1, path: 'nplanes' },
          nchannels: { type: 'child', value: 1, path: 'nchannels' },
          force_sktiff: { type: 'child', value: false, path: 'force_sktiff' },
          batch_size: { type: 'child', value: 500, path: 'batch_size' },
          do_registration: { type: 'child', value: 1, path: 'do_registration' },
        },
        isUpdated: false,
      },
    },
    flowElement: {
      flowElements: [
        {
          id: 'input_0',
          type: 'ImageFileNode',
          data: { type: 'input', label: 'hoge.tif' },
          style: { border: '1px solid #777', height: 120 },
          position: { x: 50, y: 150 },
        },
        {
          id: 'suite2p_roi_16ffq69kkc',
          type: 'AlgorithmNode',
          data: { label: 'suite2p_roi', type: 'algorithm' },
          style: { width: 180, height: 100, padding: 0, borderRadius: 0 },
          targetPosition: 'left',
          sourcePosition: 'right',
          position: { x: 657, y: 170.0246319539301 },
        },
        {
          id: 'suite2p_file_convert_m58owcejm0',
          type: 'AlgorithmNode',
          data: { label: 'suite2p_file_convert', type: 'algorithm' },
          style: { width: 180, height: 100, padding: 0, borderRadius: 0 },
          targetPosition: 'left',
          sourcePosition: 'right',
          position: { x: 413, y: 169.6170204394443 },
        },
        {
          source: 'input_0',
          sourceHandle: 'input_0--image--ImageData',
          target: 'suite2p_file_convert_m58owcejm0',
          targetHandle: 'suite2p_file_convert_m58owcejm0--image--ImageData',
          animated: false,
          style: { width: 5 },
          type: 'buttonedge',
          id: 'reactflow__edge-input_0input_0--image--ImageData-suite2p_file_convert_m58owcejm0suite2p_file_convert_m58owcejm0--image--ImageData',
        },
        {
          source: 'suite2p_file_convert_m58owcejm0',
          sourceHandle: 'suite2p_file_convert_m58owcejm0--ops--Suite2pData',
          target: 'suite2p_roi_16ffq69kkc',
          targetHandle: 'suite2p_roi_16ffq69kkc--ops--Suite2pData',
          animated: false,
          style: { width: 5 },
          type: 'buttonedge',
          id: 'reactflow__edge-suite2p_file_convert_m58owcejm0suite2p_file_convert_m58owcejm0--ops--Suite2pData-suite2p_roi_16ffq69kkcsuite2p_roi_16ffq69kkc--ops--Suite2pData',
        },
      ],
      flowPosition: { x: 0, y: 0, zoom: 0.7 },
      elementCoord: { x: 600, y: 158.6170204394443 },
    },
    inputNode: {
      input_0: {
        fileType: 'image',
        param: {},
        selectedFilePath: ['/tmp/optinist/input/hoge/hoge.tif'],
      },
    },
    nwb: {
      params: {
        session_description: {
          type: 'child',
          value: 'optinist',
          path: 'session_description',
        },
        identifier: { type: 'child', value: 'optinist', path: 'identifier' },
        experiment_description: {
          type: 'child',
          value: 'None',
          path: 'experiment_description',
        },
        device: {
          type: 'parent',
          children: {
            name: {
              type: 'child',
              value: 'Microscope device',
              path: 'device/name',
            },
            description: {
              type: 'child',
              value: 'Microscope Information',
              path: 'device/description',
            },
            manufacturer: {
              type: 'child',
              value: 'Microscope Manufacture',
              path: 'device/manufacturer',
            },
          },
        },
        optical_channel: {
          type: 'parent',
          children: {
            name: {
              type: 'child',
              value: 'OpticalChannel',
              path: 'optical_channel/name',
            },
            description: {
              type: 'child',
              value: 'optical channel',
              path: 'optical_channel/description',
            },
            emission_lambda: {
              type: 'child',
              value: 500,
              path: 'optical_channel/emission_lambda',
            },
          },
        },
        imaging_plane: {
          type: 'parent',
          children: {
            name: {
              type: 'child',
              value: 'ImagingPlane',
              path: 'imaging_plane/name',
            },
            description: {
              type: 'child',
              value: 'standard',
              path: 'imaging_plane/description',
            },
            imaging_rate: {
              type: 'child',
              value: 30,
              path: 'imaging_plane/imaging_rate',
            },
            excitation_lambda: {
              type: 'child',
              value: 600,
              path: 'imaging_plane/excitation_lambda',
            },
            indicator: {
              type: 'child',
              value: 'GCaMap',
              path: 'imaging_plane/indicator',
            },
            location: {
              type: 'child',
              value: 'V1',
              path: 'imaging_plane/location',
            },
          },
        },
        image_series: {
          type: 'parent',
          children: {
            name: {
              type: 'child',
              value: 'TwoPhotonSeries',
              path: 'image_series/name',
            },
            starting_time: {
              type: 'child',
              value: 0,
              path: 'image_series/starting_time',
            },
            starting_frame: {
              type: 'child',
              value: 0,
              path: 'image_series/starting_frame',
            },
          },
        },
        ophys: {
          type: 'parent',
          children: {
            plane_segmentation: {
              type: 'parent',
              children: {
                name: {
                  type: 'child',
                  value: 'PlaneSegmentation',
                  path: 'ophys/plane_segmentation/name',
                },
                description: {
                  type: 'child',
                  value: '',
                  path: 'ophys/plane_segmentation/description',
                },
              },
            },
          },
        },
      },
    },
    snakemake: {
      params: {
        use_conda: { type: 'child', value: true, path: 'use_conda' },
        cores: { type: 'child', value: 2, path: 'cores' },
        forceall: { type: 'child', value: false, path: 'forceall' },
        forcetargets: { type: 'child', value: true, path: 'forcetargets' },
        lock: { type: 'child', value: false, path: 'lock' },
      },
    },
  } as RootState

  const expectRunPostData = {
    // name: 'post data test', // nameはomitされているので含めない
    nwbParam: {
      session_description: {
        type: 'child',
        value: 'optinist',
        path: 'session_description',
      },
      identifier: {
        type: 'child',
        value: 'optinist',
        path: 'identifier',
      },
      experiment_description: {
        type: 'child',
        value: 'None',
        path: 'experiment_description',
      },
      device: {
        type: 'parent',
        children: {
          name: {
            type: 'child',
            value: 'Microscope device',
            path: 'device/name',
          },
          description: {
            type: 'child',
            value: 'Microscope Information',
            path: 'device/description',
          },
          manufacturer: {
            type: 'child',
            value: 'Microscope Manufacture',
            path: 'device/manufacturer',
          },
        },
      },
      optical_channel: {
        type: 'parent',
        children: {
          name: {
            type: 'child',
            value: 'OpticalChannel',
            path: 'optical_channel/name',
          },
          description: {
            type: 'child',
            value: 'optical channel',
            path: 'optical_channel/description',
          },
          emission_lambda: {
            type: 'child',
            value: 500,
            path: 'optical_channel/emission_lambda',
          },
        },
      },
      imaging_plane: {
        type: 'parent',
        children: {
          name: {
            type: 'child',
            value: 'ImagingPlane',
            path: 'imaging_plane/name',
          },
          description: {
            type: 'child',
            value: 'standard',
            path: 'imaging_plane/description',
          },
          imaging_rate: {
            type: 'child',
            value: 30,
            path: 'imaging_plane/imaging_rate',
          },
          excitation_lambda: {
            type: 'child',
            value: 600,
            path: 'imaging_plane/excitation_lambda',
          },
          indicator: {
            type: 'child',
            value: 'GCaMap',
            path: 'imaging_plane/indicator',
          },
          location: {
            type: 'child',
            value: 'V1',
            path: 'imaging_plane/location',
          },
        },
      },
      image_series: {
        type: 'parent',
        children: {
          name: {
            type: 'child',
            value: 'TwoPhotonSeries',
            path: 'image_series/name',
          },
          starting_time: {
            type: 'child',
            value: 0,
            path: 'image_series/starting_time',
          },
          starting_frame: {
            type: 'child',
            value: 0,
            path: 'image_series/starting_frame',
          },
        },
      },
      ophys: {
        type: 'parent',
        children: {
          plane_segmentation: {
            type: 'parent',
            children: {
              name: {
                type: 'child',
                value: 'PlaneSegmentation',
                path: 'ophys/plane_segmentation/name',
              },
              description: {
                type: 'child',
                value: '',
                path: 'ophys/plane_segmentation/description',
              },
            },
          },
        },
      },
    },
    snakemakeParam: {
      use_conda: { type: 'child', value: true, path: 'use_conda' },
      cores: { type: 'child', value: 2, path: 'cores' },
      forceall: { type: 'child', value: false, path: 'forceall' },
      forcetargets: { type: 'child', value: true, path: 'forcetargets' },
      lock: { type: 'child', value: false, path: 'lock' },
    },
    edgeDict: {
      'reactflow__edge-input_0input_0--image--ImageData-suite2p_file_convert_m58owcejm0suite2p_file_convert_m58owcejm0--image--ImageData':
        {
          source: 'input_0',
          sourceHandle: 'input_0--image--ImageData',
          target: 'suite2p_file_convert_m58owcejm0',
          targetHandle: 'suite2p_file_convert_m58owcejm0--image--ImageData',
          animated: false,
          style: { width: 5 },
          type: 'buttonedge',
          id: 'reactflow__edge-input_0input_0--image--ImageData-suite2p_file_convert_m58owcejm0suite2p_file_convert_m58owcejm0--image--ImageData',
        },
      'reactflow__edge-suite2p_file_convert_m58owcejm0suite2p_file_convert_m58owcejm0--ops--Suite2pData-suite2p_roi_16ffq69kkcsuite2p_roi_16ffq69kkc--ops--Suite2pData':
        {
          source: 'suite2p_file_convert_m58owcejm0',
          sourceHandle: 'suite2p_file_convert_m58owcejm0--ops--Suite2pData',
          target: 'suite2p_roi_16ffq69kkc',
          targetHandle: 'suite2p_roi_16ffq69kkc--ops--Suite2pData',
          animated: false,
          style: { width: 5 },
          type: 'buttonedge',
          id: 'reactflow__edge-suite2p_file_convert_m58owcejm0suite2p_file_convert_m58owcejm0--ops--Suite2pData-suite2p_roi_16ffq69kkcsuite2p_roi_16ffq69kkc--ops--Suite2pData',
        },
    },
    nodeDict: {
      input_0: {
        id: 'input_0',
        type: 'ImageFileNode',
        data: {
          type: 'input',
          label: 'hoge.tif',
          path: ['/tmp/optinist/input/hoge/hoge.tif'],
          param: {},
          fileType: 'image',
        },
        style: { border: '1px solid #777', height: 120 },
        position: { x: 50, y: 150 },
      },
      suite2p_roi_16ffq69kkc: {
        id: 'suite2p_roi_16ffq69kkc',
        type: 'AlgorithmNode',
        data: {
          label: 'suite2p_roi',
          type: 'algorithm',
          path: 'suite2p/suite2p_roi',
          param: {
            tau: { type: 'child', value: 1, path: 'tau' },
            fs: { type: 'child', value: 10, path: 'fs' },
            soma_crop: { type: 'child', value: true, path: 'soma_crop' },
            high_pass: { type: 'child', value: 100, path: 'high_pass' },
            sparse_mode: {
              type: 'child',
              value: true,
              path: 'sparse_mode',
            },
            max_overlap: {
              type: 'child',
              value: 0.75,
              path: 'max_overlap',
            },
            nbinned: { type: 'child', value: 5000, path: 'nbinned' },
            spatial_scale: {
              type: 'child',
              value: 0,
              path: 'spatial_scale',
            },
            threshold_scaling: {
              type: 'child',
              value: 1,
              path: 'threshold_scaling',
            },
            max_iterations: {
              type: 'child',
              value: 20,
              path: 'max_iterations',
            },
            spatial_hp_detect: {
              type: 'child',
              value: 25,
              path: 'spatial_hp_detect',
            },
            preclassify: { type: 'child', value: 0, path: 'preclassify' },
            allow_overlap: {
              type: 'child',
              value: false,
              path: 'allow_overlap',
            },
            inner_neuropil_radius: {
              type: 'child',
              value: 2,
              path: 'inner_neuropil_radius',
            },
            min_neuropil_pixels: {
              type: 'child',
              value: 350,
              path: 'min_neuropil_pixels',
            },
            neucoeff: { type: 'child', value: 0.7, path: 'neucoeff' },
          },
        },
        style: { width: 180, height: 100, padding: 0, borderRadius: 0 },
        targetPosition: 'left',
        sourcePosition: 'right',
        position: { x: 657, y: 170.0246319539301 },
      },
      suite2p_file_convert_m58owcejm0: {
        id: 'suite2p_file_convert_m58owcejm0',
        type: 'AlgorithmNode',
        data: {
          label: 'suite2p_file_convert',
          type: 'algorithm',
          path: 'suite2p/suite2p_file_convert',
          param: {
            nplanes: { type: 'child', value: 1, path: 'nplanes' },
            nchannels: { type: 'child', value: 1, path: 'nchannels' },
            force_sktiff: {
              type: 'child',
              value: false,
              path: 'force_sktiff',
            },
            batch_size: { type: 'child', value: 500, path: 'batch_size' },
            do_registration: {
              type: 'child',
              value: 1,
              path: 'do_registration',
            },
          },
        },
        style: { width: 180, height: 100, padding: 0, borderRadius: 0 },
        targetPosition: 'left',
        sourcePosition: 'right',
        position: { x: 413, y: 169.6170204394443 },
      },
    },
    forceRunList: [],
  }

  const targetStateForceRunList = {
    ...initialRootState,
    algorithmNode: {
      suite2p_file_convert_6fn2k01zph: {
        functionPath: 'suite2p/suite2p_file_convert',
        name: 'suite2p_file_convert',
        params: {
          nplanes: { type: 'child', value: 2, path: 'nplanes' },
          nchannels: { type: 'child', value: 1, path: 'nchannels' },
          force_sktiff: { type: 'child', value: true, path: 'force_sktiff' },
          batch_size: { type: 'child', value: 500, path: 'batch_size' },
          do_registration: { type: 'child', value: 1, path: 'do_registration' },
        },
        isUpdated: true,
      },
    },
    flowElement: {
      flowElements: [
        {
          id: 'input_0',
          type: 'ImageFileNode',
          data: { type: 'input', label: 'hoge.tif' },
          style: { border: '1px solid #777', height: 120 },
          position: { x: 50, y: 152 },
        },
        {
          id: 'suite2p_file_convert_6fn2k01zph',
          type: 'AlgorithmNode',
          data: { label: 'suite2p_file_convert', type: 'algorithm' },
          style: { width: 180, height: 100, padding: 0, borderRadius: 0 },
          targetPosition: 'left',
          sourcePosition: 'right',
          position: { x: 350, y: 153.52022229668373 },
        },
        {
          source: 'input_0',
          sourceHandle: 'input_0--image--ImageData',
          target: 'suite2p_file_convert_6fn2k01zph',
          targetHandle: 'suite2p_file_convert_6fn2k01zph--image--ImageData',
          animated: false,
          style: { width: 5 },
          type: 'buttonedge',
          id: 'reactflow__edge-input_0input_0--image--ImageData-suite2p_file_convert_6fn2k01zphsuite2p_file_convert_6fn2k01zph--image--ImageData',
        },
      ],
      flowPosition: { x: 0, y: 0, zoom: 0.7 },
      elementCoord: { x: 350, y: 153.52022229668373 },
    },
    inputNode: {
      input_0: {
        fileType: 'image',
        param: {},
        selectedFilePath: ['/tmp/optinist/input/hoge/hoge.tif'],
      },
    },
    nwb: { params: {} },
    snakemake: { params: {} },
  } as RootState

  const expectRunPostDataForceRunList = {
    nwbParam: {},
    snakemakeParam: {},
    edgeDict: {
      'reactflow__edge-input_0input_0--image--ImageData-suite2p_file_convert_6fn2k01zphsuite2p_file_convert_6fn2k01zph--image--ImageData':
        {
          source: 'input_0',
          sourceHandle: 'input_0--image--ImageData',
          target: 'suite2p_file_convert_6fn2k01zph',
          targetHandle: 'suite2p_file_convert_6fn2k01zph--image--ImageData',
          animated: false,
          style: { width: 5 },
          type: 'buttonedge',
          id: 'reactflow__edge-input_0input_0--image--ImageData-suite2p_file_convert_6fn2k01zphsuite2p_file_convert_6fn2k01zph--image--ImageData',
        },
    },
    nodeDict: {
      input_0: {
        id: 'input_0',
        type: 'ImageFileNode',
        data: {
          type: 'input',
          label: 'hoge.tif',
          path: ['/tmp/optinist/input/hoge/hoge.tif'],
          param: {},
          fileType: 'image',
        },
        style: { border: '1px solid #777', height: 120 },
        position: { x: 50, y: 152 },
      },
      suite2p_file_convert_6fn2k01zph: {
        id: 'suite2p_file_convert_6fn2k01zph',
        type: 'AlgorithmNode',
        data: {
          label: 'suite2p_file_convert',
          type: 'algorithm',
          path: 'suite2p/suite2p_file_convert',
          param: {
            nplanes: { type: 'child', value: 2, path: 'nplanes' },
            nchannels: { type: 'child', value: 1, path: 'nchannels' },
            force_sktiff: {
              type: 'child',
              value: true,
              path: 'force_sktiff',
            },
            batch_size: { type: 'child', value: 500, path: 'batch_size' },
            do_registration: {
              type: 'child',
              value: 1,
              path: 'do_registration',
            },
          },
        },
        style: { width: 180, height: 100, padding: 0, borderRadius: 0 },
        targetPosition: 'left',
        sourcePosition: 'right',
        position: { x: 350, y: 153.52022229668373 },
      },
    },
    forceRunList: [
      {
        nodeId: 'suite2p_file_convert_6fn2k01zph',
        name: 'suite2p_file_convert',
      },
    ],
  }
})
