import { Edge } from "reactflow"

import { expect, describe, test, beforeEach } from "@jest/globals"
import { AsyncThunk } from "@reduxjs/toolkit"

import {
  AlgorithmNodePostData,
  EdgeDict,
  InputNodePostData,
  NodeDict,
  RunPostData,
} from "api/run/Run"
import {
  run,
  runByCurrentUid,
  pollRunResult,
  cancelResult,
} from "store/slice/Pipeline/PipelineActions"
import * as selectors from "store/slice/Pipeline/PipelineSelectors"
import reducer, { initialState } from "store/slice/Pipeline/PipelineSlice"
import {
  NODE_RESULT_STATUS,
  Pipeline,
  PipelineType,
  RUN_STATUS,
} from "store/slice/Pipeline/PipelineType"
import { isStartedPipeline } from "store/slice/Pipeline/PipelineUtils"
import { RootState } from "store/store"
import { ParamMap, ParamType } from "utils/param/ParamType"

const requestId = "FmYmw6sCHA2Ll5JJfPuJN"

// Factory for reusable mock data structures
const createNodeDict = (): NodeDict => ({
  node1: {
    id: "node1",
    position: { x: 0, y: 0 },
    data: {
      type: "input",
      label: "input1",
      path: "/path/to/input",
      fileType: "image",
    } as InputNodePostData,
  },
  node2: {
    id: "node2",
    position: { x: 1, y: 1 },
    data: {
      path: "/path/to/algorithm",
      param: {
        key: {
          type: "child",
          value: "",
          path: "key",
          children: {
            key1: {
              type: "child",
              value: "value1",
              path: "key.key1",
            } as ParamType,
          },
        },
      } as ParamMap,
    } as AlgorithmNodePostData,
  },
})

const createEdgeDict = (): EdgeDict => ({
  node1: { id: "node1", type: "type1", data: { key: "value1" } } as Edge,
  node2: { id: "node2", type: "type2", data: { key: "value2" } } as Edge,
})

const createSnakemakeParam = (): ParamMap => ({
  config: {
    type: "parent",
    children: {
      cores: { type: "child", value: 4, path: "config.cores" },
      memory: { type: "child", value: "16GB", path: "config.memory" },
      restartTimes: { type: "child", value: 3, path: "config.restartTimes" },
    },
  },
})

const createNwbParam = (): ParamMap => ({
  config: {
    type: "parent",
    children: {
      filePath: {
        type: "child",
        value: "/path/to/nwb/file",
        path: "nwbParam.filePath",
      },
      options: {
        type: "child",
        value: { option1: true, option2: false },
        path: "nwbParam.options",
      },
    },
  },
})

const runPostData: RunPostData = {
  name: "record test",
  nodeDict: createNodeDict(),
  edgeDict: createEdgeDict(),
  snakemakeParam: createSnakemakeParam(),
  nwbParam: createNwbParam(),
  forceRunList: [],
}

// Helper to create a fulfilled action for a given action type and payload
const createFulfilledAction = (
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  actionType: AsyncThunk<any, any, any>,
  payload = "response data",
) => ({
  type: actionType.fulfilled.type,
  meta: {
    requestId: requestId,
    requestStatus: "fulfilled",
    arg: { runPostData },
  },
  payload,
})

const runPendingAction = {
  type: run.pending.type,
  meta: {
    requestId: requestId,
    requestStatus: "pending",
  },
}

const runRejectedAction = {
  type: run.rejected.type,
  meta: {
    requestId: requestId,
    requestStatus: "rejected",
    arg: { runPostData },
  },
  error: "error message",
}

const runByCurrentUidPendingAction = {
  type: runByCurrentUid.pending.type,
  meta: { requestId: requestId, requestStatus: "pending" },
}

const pollRunResultPendingAction = {
  type: pollRunResult.pending.type,
  meta: { requestId: requestId, requestStatus: "pending" },
}

describe("Pipeline State Test", () => {
  describe("Pipeline Run", () => {
    test(run.fulfilled.type, () => {
      const runFulfilledAction = createFulfilledAction(run)
      const targetState = reducer(
        reducer(initialState, runPendingAction),
        runFulfilledAction,
      )

      const expectState = {
        currentPipeline: { uid: "response data" },
        run: {
          runPostData,
          runResult: {},
          status: RUN_STATUS.START_SUCCESS,
          uid: "response data",
        },
        runBtn: 1,
      }
      expect(targetState).toEqual(expectState)
    })

    test(run.pending.type, () => {
      const targetState = reducer(initialState, runPendingAction)
      const expectState = {
        run: { status: RUN_STATUS.START_PENDING },
        runBtn: 1,
      }
      expect(targetState).toEqual(expectState)
    })

    test(run.rejected.type, () => {
      const targetState = reducer(
        reducer(initialState, runPendingAction),
        runRejectedAction,
      )
      const expectState = { run: { status: RUN_STATUS.START_ERROR }, runBtn: 1 }
      expect(targetState).toEqual(expectState)
    })
  })

  describe("Pipeline RunByCurrentId", () => {
    test(runByCurrentUid.fulfilled.type, () => {
      const runByCurrentUidFulfilledAction =
        createFulfilledAction(runByCurrentUid)
      const targetState = reducer(
        reducer(initialState, runByCurrentUidPendingAction),
        runByCurrentUidFulfilledAction,
      )

      const expectState = {
        currentPipeline: { uid: "response data" },
        run: {
          runPostData,
          runResult: {},
          status: RUN_STATUS.START_SUCCESS,
          uid: "response data",
        },
        runBtn: 1,
      }
      expect(targetState).toEqual(expectState)
    })

    test(runByCurrentUid.pending.type, () => {
      const targetState = reducer(initialState, runByCurrentUidPendingAction)
      const expectState = {
        run: { status: RUN_STATUS.START_PENDING },
        runBtn: 1,
      }
      expect(targetState).toEqual(expectState)
    })

    test(runByCurrentUid.rejected.type, () => {
      const runByCurrentUidRejectedAction = {
        type: runByCurrentUid.rejected.type,
        meta: {
          requestId: requestId,
          requestStatus: "rejected",
          arg: { runPostData },
        },
        error: "error message",
      }
      const targetState = reducer(
        reducer(initialState, runByCurrentUidPendingAction),
        runByCurrentUidRejectedAction,
      )
      const expectState = { run: { status: RUN_STATUS.START_ERROR }, runBtn: 1 }
      expect(targetState).toEqual(expectState)
    })
  })

  describe("Pipeline PollRunResult", () => {
    test(pollRunResult.fulfilled.type, () => {
      const pollRunResultFulfilledAction = {
        type: pollRunResult.fulfilled.type,
        meta: {
          requestId: requestId,
          requestStatus: "fulfilled",
          arg: { uid: "test-uid" },
        },
        payload: {
          node1: {
            status: "success",
            message: "Node 1 completed successfully",
            name: "Node 1",
            outputPaths: {
              output1: { path: "/path/to/output1", type: "images" },
            },
          },
          node2: { status: "failed", message: "Node 2 failed", name: "Node 2" },
        },
      }

      const mockState: Pipeline = {
        ...initialState,
        run: {
          ...initialState.run,
          status: "StartSuccess" as const, // Ensures TypeScript treats this as a literal
          runResult: {
            node1: {
              status: "success",
              name: "Node 1",
              message: "",
              outputPaths: {
                output1: { path: "./output/path", type: "image" },
              },
            },
            node2: { status: "pending", name: "Node 2", message: "" },
          },
          uid: "test-uid",
          runPostData,
        },
      }

      const targetState = reducer(
        reducer(mockState, pollRunResultPendingAction),
        pollRunResultFulfilledAction,
      )

      const expectState = {
        ...initialState,
        run: {
          ...initialState.run,
          runPostData,
          runResult: {
            node1: {
              status: "success",
              message: "Node 1 completed successfully",
              name: "Node 1",
              outputPaths: {
                output1: { path: "/path/to/output1", type: "image" },
              },
            },
            node2: {
              status: "error",
              message: "Node 2 failed",
              name: "Node 2",
            },
          },
          status: RUN_STATUS.FINISHED,
          uid: "test-uid",
        },
      }
      expect(targetState).toEqual(expectState)
    })

    test(pollRunResult.rejected.type, () => {
      const pollRunResultRejectedAction = {
        type: pollRunResult.rejected.type,
        meta: {
          requestId: requestId,
          requestStatus: "rejected",
          arg: { uid: "test-uid" },
        },
        error: "error message",
      }
      const targetState = reducer(
        reducer(initialState, pollRunResultPendingAction),
        pollRunResultRejectedAction,
      )
      const expectState = {
        ...initialState,
        run: { ...initialState.run, status: RUN_STATUS.ABORTED },
      }
      expect(targetState).toEqual(expectState)
    })
  })

  describe("Pipeline CancelResult", () => {
    const cancelResultPendingAction = {
      type: cancelResult.pending.type,
      meta: { requestId: requestId, requestStatus: "pending" },
    }

    const cancelResultRejectedAction = {
      type: cancelResult.rejected.type,
      meta: {
        requestId: requestId,
        requestStatus: "rejected",
        arg: { uid: "test-uid" },
      },
      error: "error message",
    }

    test(cancelResult.fulfilled.type, () => {
      const cancelResultFulfilledAction = {
        type: cancelResult.fulfilled.type,
        meta: {
          requestId: requestId,
          requestStatus: "fulfilled",
          arg: { uid: "test-uid" },
        },
        payload: {
          message: "Cancellation successful",
        },
      }

      const targetState = reducer(
        reducer(initialState, cancelResultPendingAction),
        cancelResultFulfilledAction,
      )

      const expectState = {
        ...initialState,
        run: {
          ...initialState.run,
          status: RUN_STATUS.CANCELED,
        },
      }
      expect(targetState).toEqual(expectState)
    })

    // Status is not changing when CancelResult is pending at PipelineSlice
    test(cancelResult.pending.type, () => {
      const targetState = reducer(initialState, cancelResultPendingAction)
      const expectState = {
        run: { status: RUN_STATUS.START_UNINITIALIZED },
        runBtn: 1,
      }
      expect(targetState).toEqual(expectState)
    })

    // Status is not changing when CancelResult is rejected at PipelineSlice
    test(cancelResult.rejected.type, () => {
      const targetState = reducer(
        reducer(initialState, cancelResultPendingAction),
        cancelResultRejectedAction,
      )
      const expectState = {
        ...initialState,
        run: { ...initialState.run, status: RUN_STATUS.START_UNINITIALIZED },
      }
      expect(targetState).toEqual(expectState)
    })
  })
})

describe("PipelineUtils isStartedPipeline", () => {
  test("returns true for status RUN_STATUS.START_SUCCESS", () => {
    const pipeline: PipelineType = {
      status: RUN_STATUS.START_SUCCESS,
      uid: "test-uid",
      runPostData: runPostData,
      runResult: {},
    }
    expect(isStartedPipeline(pipeline)).toBe(true)
  })

  test("returns true for status RUN_STATUS.FINISHED", () => {
    const pipeline: PipelineType = {
      status: RUN_STATUS.FINISHED,
      uid: "test-uid",
      runPostData: runPostData,
      runResult: {},
    }
    expect(isStartedPipeline(pipeline)).toBe(true)
  })

  test("returns true for status RUN_STATUS.ABORTED", () => {
    const pipeline: PipelineType = {
      status: RUN_STATUS.ABORTED,
      uid: "test-uid",
      runPostData: runPostData,
      runResult: {},
    }
    expect(isStartedPipeline(pipeline)).toBe(true)
  })

  test("returns false for status other than START_SUCCESS, FINISHED, or ABORTED", () => {
    const pipeline: PipelineType = {
      status: RUN_STATUS.START_PENDING, // Example of a non-started status
      // add other necessary fields for PipelineType if needed
    }
    expect(isStartedPipeline(pipeline)).toBe(false)
  })
})

describe("Pipeline Selectors", () => {
  let initialState: RootState

  beforeEach(() => {
    initialState = {
      pipeline: {
        currentPipeline: { uid: "123" },
        run: {
          runResult: {
            node1: {
              status: "pending",
              message: "Processing...",
              nodeResult: {
                status: NODE_RESULT_STATUS.PENDING,
              },
            },
            node2: {
              status: "success",
              message: "Completed",
              outputPaths: { key1: { path: "/output/file1", type: "csv" } },
            },
          },
          status: RUN_STATUS.START_SUCCESS,
        },
        runBtn: false,
      },
      experiments: {
        status: "fulfilled",
        experimentList: {
          "123": { name: "Pipeline 123" },
        },
      },
    } as unknown as RootState
  })

  test("selectPipelineLatestUid should return the latest pipeline UID", () => {
    const result = selectors.selectPipelineLatestUid(initialState)
    expect(result).toBe("123")
  })

  test("selectCurrentPipelineName should return the current pipeline name", () => {
    const result = selectors.selectCurrentPipelineName(initialState)
    expect(result).toBe("Pipeline 123")
  })

  test("selectRunResultPendingList should return pending nodes", () => {
    const result = selectors.selectRunResultPendingList(initialState)
    const expectedValue = [
      {
        message: "Processing...",
        nodeResult: { status: "pending" },
        status: "pending",
      },
    ]
    expect(result).toEqual(expectedValue)
  })

  test("selectRunResultPendingNodeIdList should return pending node IDs", () => {
    const result = selectors.selectRunResultPendingNodeIdList(initialState)
    expect(result).toEqual(["node1"])
  })

  test("selectPipelineIsStartedSuccess should return true if pipeline is started successfully", () => {
    const result = selectors.selectPipelineIsStartedSuccess(initialState)
    expect(result).toBe(true)
  })

  test("selectPipelineNodeResultOutputFilePath should return the output file path for a given node and key", () => {
    const result = selectors.selectPipelineNodeResultOutputFilePath(
      "node2",
      "key1",
    )(initialState)
    expect(result).toBe("/output/file1")
  })

  test("selectPipelineNodeResultOutputFileDataType should return the data type for a given node and key", () => {
    const result = selectors.selectPipelineNodeResultOutputFileDataType(
      "node2",
      "key1",
    )(initialState)
    expect(result).toBe("csv")
  })

  test("selectPipelineNodeResultMessage should return the message for a given node", () => {
    const result =
      selectors.selectPipelineNodeResultMessage("node1")(initialState)
    expect(result).toBe("Processing...")
  })

  test("selectPipelineNodeResultOutputKeyList should return the output keys for a given node", () => {
    const result =
      selectors.selectPipelineNodeResultOutputKeyList("node2")(initialState)
    expect(result).toEqual(["key1"])
  })

  test("selectPipelineNodeResultSuccessList should return nodes with success status", () => {
    const result = selectors.selectPipelineNodeResultSuccessList(initialState)
    expect(result).toEqual([
      {
        nodeId: "node2",
        nodeResult: {
          status: "success",
          message: "Completed",
          outputPaths: { key1: { path: "/output/file1", type: "csv" } },
        },
      },
    ])
  })
})
