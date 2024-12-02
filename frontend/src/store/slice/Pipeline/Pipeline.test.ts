import { expect, describe, test, beforeEach } from "@jest/globals"
import { AnyAction } from "@reduxjs/toolkit"

import {
  run,
  runByCurrentUid,
  pollRunResult,
  cancelResult,
} from "store/slice/Pipeline/PipelineActions"
import * as selectors from "store/slice/Pipeline/PipelineSelectors"
import reducer, { initialState } from "store/slice/Pipeline/PipelineSlice"
import {
  createFulfilledWithRunPostDataAction,
  createPendingAction,
  createPendingWithRequestIdAction,
  createRejectedAction,
  createRunResultAndCacelFulfilledAction,
  createRunResultAndCancelRejectedAction,
  runPostData,
} from "store/slice/Pipeline/PipelineTestUtils"
import {
  NODE_RESULT_STATUS,
  Pipeline,
  PipelineType,
  RUN_STATUS,
} from "store/slice/Pipeline/PipelineType"
import { isStartedPipeline } from "store/slice/Pipeline/PipelineUtils"
import { RootState } from "store/store"

const pollRunResultPayload = {
  node1: {
    status: "success",
    message: "Node 1 completed successfully",
    name: "Node 1",
    outputPaths: {
      output1: { path: "/path/to/output1", type: "images" },
    },
  },
  node2: { status: "failed", message: "Node 2 failed", name: "Node 2" },
}

describe("Pipeline State Test", () => {
  describe("Pipeline Run", () => {
    test(run.fulfilled.type, () => {
      const runFulfilledAction = createFulfilledWithRunPostDataAction(
        run,
      ) as AnyAction
      const targetState = reducer(
        reducer(initialState, createPendingWithRequestIdAction(run)),
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
      const targetState = reducer(
        initialState,
        createPendingWithRequestIdAction(run),
      )
      const expectState = {
        run: { status: RUN_STATUS.START_PENDING },
        runBtn: 1,
      }
      expect(targetState).toEqual(expectState)
    })

    test(run.rejected.type, () => {
      const targetState = reducer(
        reducer(initialState, createPendingWithRequestIdAction(run)),
        createRejectedAction(run),
      )
      const expectState = { run: { status: RUN_STATUS.START_ERROR }, runBtn: 1 }
      expect(targetState).toEqual(expectState)
    })
  })

  describe("Pipeline RunByCurrentId", () => {
    test(runByCurrentUid.fulfilled.type, () => {
      const runByCurrentUidFulfilledAction =
        createFulfilledWithRunPostDataAction(runByCurrentUid)
      const targetState = reducer(
        reducer(initialState, createPendingAction(runByCurrentUid)),
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
      const targetState = reducer(
        initialState,
        createPendingAction(runByCurrentUid),
      )
      const expectState = {
        run: { status: RUN_STATUS.START_PENDING },
        runBtn: 1,
      }
      expect(targetState).toEqual(expectState)
    })

    test(runByCurrentUid.rejected.type, () => {
      const targetState = reducer(
        reducer(initialState, createPendingAction(runByCurrentUid)),
        createRejectedAction(runByCurrentUid),
      )
      const expectState = { run: { status: RUN_STATUS.START_ERROR }, runBtn: 1 }
      expect(targetState).toEqual(expectState)
    })
  })

  describe("Pipeline PollRunResult", () => {
    test(pollRunResult.fulfilled.type, () => {
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
        reducer(mockState, createPendingAction(pollRunResult)),
        createRunResultAndCacelFulfilledAction(
          pollRunResult,
          pollRunResultPayload,
        ),
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
      const targetState = reducer(
        reducer(initialState, createPendingAction(pollRunResult)),
        createRunResultAndCancelRejectedAction(pollRunResult, "error message"),
      )
      const expectState = {
        ...initialState,
        run: { ...initialState.run, status: RUN_STATUS.ABORTED },
      }
      expect(targetState).toEqual(expectState)
    })
  })

  describe("Pipeline CancelResult", () => {
    test(cancelResult.fulfilled.type, () => {
      const targetState = reducer(
        reducer(initialState, createPendingAction(cancelResult)),
        createRunResultAndCacelFulfilledAction(
          cancelResult,
          "Cancellation successful",
        ),
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
      const targetState = reducer(
        initialState,
        createPendingWithRequestIdAction(cancelResult),
      )
      const expectState = {
        run: { status: RUN_STATUS.START_UNINITIALIZED },
        runBtn: 1,
      }
      expect(targetState).toEqual(expectState)
    })

    // Status is not changing when CancelResult is rejected at PipelineSlice
    test(cancelResult.rejected.type, () => {
      const targetState = reducer(
        reducer(initialState, createPendingWithRequestIdAction(cancelResult)),
        createRunResultAndCancelRejectedAction(cancelResult, "error message"),
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
