import { expect, describe, test } from "@jest/globals"

import { RUN_BTN_OPTIONS } from "store/slice/Pipeline/PipelineType"
import { reproduceWorkflow } from "store/slice/Workflow/WorkflowActions"
import { store, rootReducer } from "store/store"

describe("RecordReproduce", () => {
  const initialState = store.getState()

  const reproduceWorkflowPayload = {
    workspace_id: 1,
    unique_id: "96844a59",
    name: "dummy",
    started_at: "2023-10-13 15:03:13",
    finished_at: "2023-10-13 15:03:59",
    success: "success",
    hasNWB: true,
    function: {
      input_0: {
        unique_id: "input_0",
        name: "hoge.tif",
        success: "success",
        hasNWB: false,
        message: null,
        outputPaths: null,
        started_at: "2023-,10-13 15:03:13",
        finished_at: "2023-10-13 15:03:13",
      },
      dummy_image2image_c8tqfxw0mq: {
        unique_id: "dummy_image2image_c8tqfxw0mq",
        name: "dummy_image2image",
        success: "success",
        hasNWB: false,
        message: "success",
        outputPaths: {
          dummyImage: {
            path: "/tmp/optinist/output/96844a59/dummy_image2image_c8tqfxw0mq/image.tif",
            type: "images",
            max_index: 1,
          },
        },
        started_at: "2023-10-13 15:03:14",
        finished_at: "2023-10-13 15:03:40",
      },
      dummy_image2image8time_4mrz8h7hyk: {
        unique_id: "dummy_image2image8time_4mrz8h7hyk",
        name: "dummy_image2image8time",
        success: "success",
        hasNWB: false,
        outputPaths: {
          dummyImage2: {
            path: "/tmp/optinist/output/96844a59/dummy_image2image8time_4mrz8h7hyk/image2.tif",
            type: "images",
            max_index: 1,
          },
        },
        started_at: "2023-10-13 15:03:40",
        finished_at: "2023-10-13 15:03:59",
      },
    },
    nodeDict: {
      input_0: {
        id: "input_0",
        type: "ImageFileNode",
        data: {
          label: "hoge.tif",
          param: {},
          path: ["/tmp/optinist/input/hoge/hoge.tif"],
          type: "input",
          fileType: "image",
          hdf5Path: null,
        },
        position: { x: 51, y: 150 },
        style: {
          border: "1px solid #777",
          height: 120,
          padding: null,
          width: null,
          borderRadius: null,
        },
      },
      dummy_image2image_c8tqfxw0mq: {
        id: "dummy_image2image_c8tqfxw0mq",
        type: "AlgorithmNode",
        data: {
          label: "dummy_image2image",
          param: { sample: { path: "sample", type: "child", value: "test" } },
          path: "dummy/dummy_image2image",
          type: "algorithm",
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
      dummy_image2image8time_4mrz8h7hyk: {
        id: "dummy_image2image8time_4mrz8h7hyk",
        type: "AlgorithmNode",
        data: {
          label: "dummy_image2image8time",
          param: {},
          path: "dummy/dummy_image2image8time",
          type: "algorithm",
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
    },
    edgeDict: {
      "reactflow__edge-dummy_image2image_c8tqfxw0mqdummy_image2image_c8tqfxw0mq--image2image--ImageData-dummy_image2image8time_4mrz8h7hykdummy_image2image8time_4mrz8h7hyk--image1--ImageData":
        {
          id: "reactflow__edge-dummy_image2image_c8tqfxw0mqdummy_image2image_c8tqfxw0mq--image2image--ImageData-dummy_image2image8time_4mrz8h7hykdummy_image2image8time_4mrz8h7hyk--image1--ImageData",
          type: "buttonedge",
          animated: false,
          source: "dummy_image2image_c8tqfxw0mq",
          sourceHandle: "dummy_image2image_c8tqfxw0mq--image2image--ImageData",
          target: "dummy_image2image8time_4mrz8h7hyk",
          targetHandle: "dummy_image2image8time_4mrz8h7hyk--image1--ImageData",
          style: {
            border: null,
            height: null,
            padding: null,
            width: 5,
            borderRadius: null,
          },
        },
      "reactflow__edge-input_0input_0--image--ImageData-dummy_image2image_c8tqfxw0mqdummy_image2image_c8tqfxw0mq--image--ImageData":
        {
          id: "reactflow__edge-input_0input_0--image--ImageData-dummy_image2image_c8tqfxw0mqdummy_image2image_c8tqfxw0mq--image--ImageData",
          type: "buttonedge",
          animated: false,
          source: "input_0",
          sourceHandle: "input_0--image--ImageData",
          target: "dummy_image2image_c8tqfxw0mq",
          targetHandle: "dummy_image2image_c8tqfxw0mq--image--ImageData",
          style: {
            border: null,
            height: null,
            padding: null,
            width: 5,
            borderRadius: null,
          },
        },
    },
  }

  const uid = "96844a59"

  const expectState = {
    ...initialState,
    algorithmNode: {
      dummy_image2image_c8tqfxw0mq: {
        name: "dummy_image2image",
        functionPath: "dummy/dummy_image2image",
        params: { sample: { path: "sample", type: "child", value: "test" } },
        originalValue: false,
      },
      dummy_image2image8time_4mrz8h7hyk: {
        name: "dummy_image2image8time",
        functionPath: "dummy/dummy_image2image8time",
        params: {},
        originalValue: false,
      },
    },
    flowElement: {
      flowEdges: [
        {
          id: "reactflow__edge-dummy_image2image_c8tqfxw0mqdummy_image2image_c8tqfxw0mq--image2image--ImageData-dummy_image2image8time_4mrz8h7hykdummy_image2image8time_4mrz8h7hyk--image1--ImageData",
          type: "buttonedge",
          animated: false,
          source: "dummy_image2image_c8tqfxw0mq",
          sourceHandle: "dummy_image2image_c8tqfxw0mq--image2image--ImageData",
          target: "dummy_image2image8time_4mrz8h7hyk",
          targetHandle: "dummy_image2image8time_4mrz8h7hyk--image1--ImageData",
          style: {
            border: null,
            height: null,
            padding: null,
            width: 5,
            borderRadius: null,
          },
        },
        {
          id: "reactflow__edge-input_0input_0--image--ImageData-dummy_image2image_c8tqfxw0mqdummy_image2image_c8tqfxw0mq--image--ImageData",
          type: "buttonedge",
          animated: false,
          source: "input_0",
          sourceHandle: "input_0--image--ImageData",
          target: "dummy_image2image_c8tqfxw0mq",
          targetHandle: "dummy_image2image_c8tqfxw0mq--image--ImageData",
          style: {
            border: null,
            height: null,
            padding: null,
            width: 5,
            borderRadius: null,
          },
        },
      ],
      flowNodes: [
        {
          data: { label: "hoge.tif", type: "input" },
          id: "input_0",
          type: "ImageFileNode",
          position: { x: 51, y: 150 },
          style: {
            border: "1px solid #777",
            height: 140,
            width: 250,
          },
        },
        {
          id: "dummy_image2image_c8tqfxw0mq",
          type: "AlgorithmNode",
          data: { label: "dummy_image2image", type: "algorithm" },
          position: { x: 350, y: 151.3534781075913 },
          style: {
            border: "1px solid #777",
            height: 140,
            width: 250,
            padding: 0,
            borderRadius: 0,
          },
        },
        {
          id: "dummy_image2image8time_4mrz8h7hyk",
          type: "AlgorithmNode",
          data: { label: "dummy_image2image8time", type: "algorithm" },
          position: { x: 600, y: 164.03341976235507 },
          style: {
            border: "1px solid #777",
            height: 140,
            width: 250,
            padding: 0,
            borderRadius: 0,
          },
        },
      ],
      flowPosition: [0, 0, 0.7],
      elementCoord: { x: 100, y: 150 },
    },
    inputNode: {
      input_0: {
        fileType: "image",
        selectedFilePath: ["/tmp/optinist/input/hoge/hoge.tif"],
        param: {},
      },
    },
    nwb: { params: {} },
    snakemake: { params: {} },
    pipeline: {
      run: {
        status: "Finished",
        runPostData: {
          edgeDict: {
            "reactflow__edge-dummy_image2image_c8tqfxw0mqdummy_image2image_c8tqfxw0mq--image2image--ImageData-dummy_image2image8time_4mrz8h7hykdummy_image2image8time_4mrz8h7hyk--image1--ImageData":
              {
                animated: false,
                id: "reactflow__edge-dummy_image2image_c8tqfxw0mqdummy_image2image_c8tqfxw0mq--image2image--ImageData-dummy_image2image8time_4mrz8h7hykdummy_image2image8time_4mrz8h7hyk--image1--ImageData",
                source: "dummy_image2image_c8tqfxw0mq",
                sourceHandle:
                  "dummy_image2image_c8tqfxw0mq--image2image--ImageData",
                style: {
                  border: null,
                  borderRadius: null,
                  height: null,
                  padding: null,
                  width: 5,
                },
                target: "dummy_image2image8time_4mrz8h7hyk",
                targetHandle:
                  "dummy_image2image8time_4mrz8h7hyk--image1--ImageData",
                type: "buttonedge",
              },
            "reactflow__edge-input_0input_0--image--ImageData-dummy_image2image_c8tqfxw0mqdummy_image2image_c8tqfxw0mq--image--ImageData":
              {
                animated: false,
                id: "reactflow__edge-input_0input_0--image--ImageData-dummy_image2image_c8tqfxw0mqdummy_image2image_c8tqfxw0mq--image--ImageData",
                source: "input_0",
                sourceHandle: "input_0--image--ImageData",
                style: {
                  border: null,
                  borderRadius: null,
                  height: null,
                  padding: null,
                  width: 5,
                },
                target: "dummy_image2image_c8tqfxw0mq",
                targetHandle: "dummy_image2image_c8tqfxw0mq--image--ImageData",
                type: "buttonedge",
              },
          },
          forceRunList: [],
          name: "dummy",
          nodeDict: {
            dummy_image2image8time_4mrz8h7hyk: {
              data: {
                fileType: null,
                hdf5Path: null,
                label: "dummy_image2image8time",
                param: {},
                path: "dummy/dummy_image2image8time",
                type: "algorithm",
              },
              id: "dummy_image2image8time_4mrz8h7hyk",
              position: {
                x: 600,
                y: 164.03341976235507,
              },
              style: {
                border: null,
                borderRadius: 0,
                height: 100,
                padding: 0,
                width: 180,
              },
              type: "AlgorithmNode",
            },
            dummy_image2image_c8tqfxw0mq: {
              data: {
                fileType: null,
                hdf5Path: null,
                label: "dummy_image2image",
                param: {
                  sample: {
                    path: "sample",
                    type: "child",
                    value: "test",
                  },
                },
                path: "dummy/dummy_image2image",
                type: "algorithm",
              },
              id: "dummy_image2image_c8tqfxw0mq",
              position: {
                x: 350,
                y: 151.3534781075913,
              },
              style: {
                border: null,
                borderRadius: 0,
                height: 100,
                padding: 0,
                width: 180,
              },
              type: "AlgorithmNode",
            },
            input_0: {
              data: {
                fileType: "image",
                hdf5Path: null,
                label: "hoge.tif",
                param: {},
                path: ["/tmp/optinist/input/hoge/hoge.tif"],
                type: "input",
              },
              id: "input_0",
              position: {
                x: 51,
                y: 150,
              },
              style: {
                border: "1px solid #777",
                borderRadius: null,
                height: 120,
                padding: null,
                width: null,
              },
              type: "ImageFileNode",
            },
          },
          nwbParam: {},
          snakemakeParam: {},
        },
        runResult: {
          dummy_image2image8time_4mrz8h7hyk: {
            message: "",
            name: "dummy_image2image8time",
            outputPaths: {
              dummyImage2: {
                path: "/tmp/optinist/output/96844a59/dummy_image2image8time_4mrz8h7hyk/image2.tif",
                type: "image",
              },
            },
            status: "success",
          },
          dummy_image2image_c8tqfxw0mq: {
            message: "success",
            name: "dummy_image2image",
            outputPaths: {
              dummyImage: {
                path: "/tmp/optinist/output/96844a59/dummy_image2image_c8tqfxw0mq/image.tif",
                type: "image",
              },
            },
            status: "success",
          },
          input_0: {
            message: "",
            name: "hoge.tif",
            status: "error",
          },
        },
        uid: "96844a59",
      },
      runBtn: RUN_BTN_OPTIONS.RUN_ALREADY,
      currentPipeline: { uid: uid },
    },
    workspace: {
      listUserShare: undefined,
      loading: false,
      workspace: { items: [], limit: 50, offset: 0, total: 0 },
      currentWorkspace: { selectedTab: 0, workspaceId: 1 },
    },
  }

  test(reproduceWorkflow.fulfilled.type, () => {
    const targetState = rootReducer(initialState, {
      type: reproduceWorkflow.fulfilled.type,
      payload: reproduceWorkflowPayload,
      meta: {
        arg: { workspaceId: 1, uid: uid },
        requestId: "9rVf2XzPaoBIbRtyVl-Zk",
        requestStatus: "fulfilled",
      },
    })
    expect(targetState).toEqual(expectState)
  })
})
