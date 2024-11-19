import { expect, describe, test } from "@jest/globals"

import { run } from "store/slice/Pipeline/PipelineActions"
import reducer, { initialState } from "store/slice/Pipeline/PipelineSlice"

describe("Pipeline", () => {
  const runPendingAction = {
    type: run.pending.type,
    meta: {
      requestId: "FmYmw6sCHA2Ll5JJfPuJN",
      requestStatus: "pending",
    },
  }
  const runFulfilledAction = {
    type: run.fulfilled.type,
    meta: {
      requestId: "FmYmw6sCHA2Ll5JJfPuJN",
      requestStatus: "pending",
      arg: {
        runPostData: {
          name: "record test",
          nodeDict: {
            node1: {
              id: "node1",
              type: "type1",
              data: {
                key: "value1",
              },
            },
            node2: {
              id: "node2",
              type: "type2",
              data: {
                key: "value2",
              },
            },
          },
          edgeDict: {
            node1: {
              id: "node1",
              type: "type1",
              data: {
                key: "value1",
              },
            },
            node2: {
              id: "node2",
              type: "type2",
              data: {
                key: "value2",
              },
            },
          },
          snakemakeParam: {
            type: "parent",
            children: {
              config: {
                type: "parent",
                children: {
                  cores: {
                    type: "child",
                    value: 4,
                    path: "config.cores",
                  },
                  memory: {
                    type: "child",
                    value: "16GB",
                    path: "config.memory",
                  },
                  restartTimes: {
                    type: "child",
                    value: 3,
                    path: "config.restartTimes",
                  },
                },
              },
            },
          },
          nwbParam: {
            type: "parent",
            children: {
              filePath: {
                type: "child",
                value: "/path/to/nwb/file",
                path: "nwbParam.filePath",
              },
              options: {
                type: "child",
                value: {
                  option1: true,
                  option2: false,
                },
                path: "nwbParam.options",
              },
            },
          },
          forceRunList: [],
        },
      },
    },
    payload: "response data",
  }

  test(run.fulfilled.type, () => {
    const state = reducer(
      reducer(initialState, runPendingAction),
      runFulfilledAction,
    )
    const expectState = {
      currentPipeline: {
        uid: "response data",
      },
      run: {
        runPostData: {
          edgeDict: {
            node1: {
              id: "node1",
              type: "type1",
              data: {
                key: "value1",
              },
            },
            node2: {
              id: "node2",
              type: "type2",
              data: {
                key: "value2",
              },
            },
          },
          forceRunList: [],
          name: "record test",
          nwbParam: {
            type: "parent",
            children: {
              filePath: {
                type: "child",
                value: "/path/to/nwb/file",
                path: "nwbParam.filePath",
              },
              options: {
                type: "child",
                value: {
                  option1: true,
                  option2: false,
                },
                path: "nwbParam.options",
              },
            },
          },
          snakemakeParam: {
            type: "parent",
            children: {
              config: {
                type: "parent",
                children: {
                  cores: {
                    type: "child",
                    value: 4,
                    path: "config.cores",
                  },
                  memory: {
                    type: "child",
                    value: "16GB",
                    path: "config.memory",
                  },
                  restartTimes: {
                    type: "child",
                    value: 3,
                    path: "config.restartTimes",
                  },
                },
              },
            },
          },
          nodeDict: {
            node1: {
              id: "node1",
              type: "type1",
              data: {
                key: "value1",
              },
            },
            node2: {
              id: "node2",
              type: "type2",
              data: {
                key: "value2",
              },
            },
          },
        },
        runResult: {},
        status: "StartSuccess",
        uid: "response data",
      },
      runBtn: 1,
    }
    expect(state).toEqual(expectState)
  })
})
