import { store, rootReducer } from 'store/store'
import { importExperimentByUid } from 'store/slice/Experiments/ExperimentsActions'
import { RUN_BTN_OPTIONS } from 'store/slice/Pipeline/PipelineType'

describe('RecordReproduce', () => {
  const initialState = store.getState()

  const importExperimentByUidPayload = {
    nodeList: [
      {
        id: 'input_0',
        type: 'ImageFileNode',
        data: {
          label: 'hoge.tif',
          param: {},
          path: ['/tmp/optinist/input/hoge/hoge.tif'],
          type: 'input',
          fileType: 'image',
          hdf5Path: null,
        },
        position: { x: 51, y: 150 },
        style: {
          border: '1px solid #777',
          height: 120,
          padding: null,
          width: null,
          borderRadius: null,
        },
      },
      {
        id: 'dummy_image2image_c8tqfxw0mq',
        type: 'AlgorithmNode',
        data: {
          label: 'dummy_image2image',
          param: { sample: { path: 'sample', type: 'child', value: 'test' } },
          path: 'dummy/dummy_image2image',
          type: 'algorithm',
          fileType: null,
          hdf5Path: null,
        },
        position: { x: 350, y: 151.3534781075913 },
        style: {
          border: null,
          height: 100,
          padding: 0,
          width: 180,
          borderRadius: 0,
        },
      },
      {
        id: 'dummy_image2image8time_4mrz8h7hyk',
        type: 'AlgorithmNode',
        data: {
          label: 'dummy_image2image8time',
          param: {},
          path: 'dummy/dummy_image2image8time',
          type: 'algorithm',
          fileType: null,
          hdf5Path: null,
        },
        position: { x: 600, y: 164.03341976235507 },
        style: {
          border: null,
          height: 100,
          padding: 0,
          width: 180,
          borderRadius: 0,
        },
      },
    ],
    edgeList: [
      {
        id: 'reactflow__edge-dummy_image2image_c8tqfxw0mqdummy_image2image_c8tqfxw0mq--image2image--ImageData-dummy_image2image8time_4mrz8h7hykdummy_image2image8time_4mrz8h7hyk--image1--ImageData',
        type: 'buttonedge',
        animated: false,
        source: 'dummy_image2image_c8tqfxw0mq',
        sourceHandle: 'dummy_image2image_c8tqfxw0mq--image2image--ImageData',
        target: 'dummy_image2image8time_4mrz8h7hyk',
        targetHandle: 'dummy_image2image8time_4mrz8h7hyk--image1--ImageData',
        style: {
          border: null,
          height: null,
          padding: null,
          width: 5,
          borderRadius: null,
        },
      },
      {
        id: 'reactflow__edge-input_0input_0--image--ImageData-dummy_image2image_c8tqfxw0mqdummy_image2image_c8tqfxw0mq--image--ImageData',
        type: 'buttonedge',
        animated: false,
        source: 'input_0',
        sourceHandle: 'input_0--image--ImageData',
        target: 'dummy_image2image_c8tqfxw0mq',
        targetHandle: 'dummy_image2image_c8tqfxw0mq--image--ImageData',
        style: {
          border: null,
          height: null,
          padding: null,
          width: 5,
          borderRadius: null,
        },
      },
    ],
  }

  const uid = '96844a59'

  const expectState = {
    ...initialState,
    algorithmNode: {
      dummy_image2image_c8tqfxw0mq: {
        name: 'dummy_image2image',
        functionPath: 'dummy/dummy_image2image',
        params: { sample: { path: 'sample', type: 'child', value: 'test' } },
        isUpdated: false,
      },
      dummy_image2image8time_4mrz8h7hyk: {
        name: 'dummy_image2image8time',
        functionPath: 'dummy/dummy_image2image8time',
        params: {},
        isUpdated: false,
      },
    },
    flowElement: {
      flowElements: [
        {
          id: 'input_0',
          type: 'ImageFileNode',
          data: { label: 'hoge.tif', type: 'input' },
          position: { x: 51, y: 150 },
          style: {
            border: '1px solid #777',
            height: 120,
            padding: null,
            width: null,
            borderRadius: null,
          },
        },
        {
          id: 'dummy_image2image_c8tqfxw0mq',
          type: 'AlgorithmNode',
          data: { label: 'dummy_image2image', type: 'algorithm' },
          position: { x: 350, y: 151.3534781075913 },
          style: {
            border: null,
            height: 100,
            padding: 0,
            width: 180,
            borderRadius: 0,
          },
        },
        {
          id: 'dummy_image2image8time_4mrz8h7hyk',
          type: 'AlgorithmNode',
          data: { label: 'dummy_image2image8time', type: 'algorithm' },
          position: { x: 600, y: 164.03341976235507 },
          style: {
            border: null,
            height: 100,
            padding: 0,
            width: 180,
            borderRadius: 0,
          },
        },
        {
          id: 'reactflow__edge-dummy_image2image_c8tqfxw0mqdummy_image2image_c8tqfxw0mq--image2image--ImageData-dummy_image2image8time_4mrz8h7hykdummy_image2image8time_4mrz8h7hyk--image1--ImageData',
          type: 'buttonedge',
          animated: false,
          source: 'dummy_image2image_c8tqfxw0mq',
          sourceHandle: 'dummy_image2image_c8tqfxw0mq--image2image--ImageData',
          target: 'dummy_image2image8time_4mrz8h7hyk',
          targetHandle: 'dummy_image2image8time_4mrz8h7hyk--image1--ImageData',
          style: {
            border: null,
            height: null,
            padding: null,
            width: 5,
            borderRadius: null,
          },
        },
        {
          id: 'reactflow__edge-input_0input_0--image--ImageData-dummy_image2image_c8tqfxw0mqdummy_image2image_c8tqfxw0mq--image--ImageData',
          type: 'buttonedge',
          animated: false,
          source: 'input_0',
          sourceHandle: 'input_0--image--ImageData',
          target: 'dummy_image2image_c8tqfxw0mq',
          targetHandle: 'dummy_image2image_c8tqfxw0mq--image--ImageData',
          style: {
            border: null,
            height: null,
            padding: null,
            width: 5,
            borderRadius: null,
          },
        },
      ],
      flowPosition: { x: 0, y: 0, zoom: 0.7 },
      elementCoord: { x: 100, y: 150 },
    },
    inputNode: {
      input_0: {
        fileType: 'image',
        selectedFilePath: ['/tmp/optinist/input/hoge/hoge.tif'],
        param: {},
      },
    },
    nwb: { params: {} },
    snakemake: { params: {} },
    pipeline: {
      run: { status: 'StartUninitialized' },
      runBtn: RUN_BTN_OPTIONS.RUN_ALREADY,
      currentPipeline: { uid: uid },
    },
  }

  test(importExperimentByUid.fulfilled.type, () => {
    const targetState = rootReducer(initialState, {
      type: importExperimentByUid.fulfilled.type,
      payload: importExperimentByUidPayload,
      meta: {
        arg: uid,
        requestId: '9rVf2XzPaoBIbRtyVl-Zk',
        requestStatus: 'fulfilled',
      },
    })
    expect(targetState).toEqual(expectState)
  })
})
