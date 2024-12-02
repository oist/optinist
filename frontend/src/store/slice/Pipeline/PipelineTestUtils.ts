/* eslint-disable @typescript-eslint/no-explicit-any */
/* eslint-disable @typescript-eslint/no-unused-vars */
import { Edge } from "reactflow"

import { AsyncThunk, AnyAction } from "@reduxjs/toolkit"

import {
  NodeDict,
  InputNodePostData,
  AlgorithmNodePostData,
  EdgeDict,
  RunPostData,
} from "api/run/Run"
import { ParamType, ParamMap } from "utils/param/ParamType"

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

export const runPostData: RunPostData = {
  name: "record test",
  nodeDict: createNodeDict(),
  edgeDict: createEdgeDict(),
  snakemakeParam: createSnakemakeParam(),
  nwbParam: createNwbParam(),
  forceRunList: [],
}

export const pollRunResultPayload = {
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

export const createFulfilledWithRunPostDataAction = (
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

export const createPendingWithRequestIdAction = (
  actionType: AsyncThunk<any, any, any>,
): AnyAction => ({
  type: actionType.pending.type,
  meta: {
    requestId: requestId,
    requestStatus: "pending",
  },
})

export const createRejectedAction = (
  actionType: AsyncThunk<any, any, any>,
  error = "error message",
) => ({
  type: actionType.rejected.type,
  meta: {
    requestId: requestId,
    requestStatus: "rejected",
    arg: { runPostData },
  },
  error,
})

export const createPendingAction = (
  actionType: AsyncThunk<any, any, any>,
): AnyAction => ({
  type: actionType.pending.type,
  meta: {
    requestStatus: "pending",
  },
})

export const createRunResultAndCancelFulfilledAction = (
  actionType: AsyncThunk<any, any, any>,
  payload: any,
) => ({
  type: actionType.fulfilled.type,
  meta: {
    requestId: requestId,
    requestStatus: "fulfilled",
    arg: { uid: "test-uid" },
  },
  payload,
})

export const createRunResultAndCancelRejectedAction = (
  actionType: AsyncThunk<any, any, any>,
  errorMessage: string,
) => ({
  type: actionType.rejected.type,
  meta: {
    requestId: requestId,
    requestStatus: "rejected",
    arg: { uid: "test-uid" },
  },
  error: errorMessage,
})
