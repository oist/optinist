import { createSlice, isAnyOf, PayloadAction } from "@reduxjs/toolkit"

import { isInputNodePostData } from "api/run/RunUtils"
import { INITIAL_IMAGE_ELEMENT_ID } from "const/flowchart"
import { WORKSPACE_TYPE } from "const/Workspace"
import { FileNodeFactory } from "factories/FileNodeFactory"
import { uploadFile } from "store/slice/FileUploader/FileUploaderActions"
import { addInputNode } from "store/slice/FlowElement/FlowElementActions"
import {
  clearFlowElements,
  deleteFlowNodes,
  deleteFlowNodeById,
} from "store/slice/FlowElement/FlowElementSlice"
import { NODE_TYPE_SET } from "store/slice/FlowElement/FlowElementType"
import { setInputNodeFilePath } from "store/slice/InputNode/InputNodeActions"
import {
  CsvInputParamType,
  FILE_TYPE_SET,
  HDF5InputNode,
  InputNode,
  InputNodeType,
  INPUT_NODE_SLICE_NAME,
  MatlabInputNode,
  WORKSPACE_TYPE_KEY,
} from "store/slice/InputNode/InputNodeType"
import {
  isCsvInputNode,
  isHDF5InputNode,
  isMatlabInputNode,
} from "store/slice/InputNode/InputNodeUtils"
import {
  reproduceWorkflow,
  importWorkflowConfig,
  fetchWorkflow,
} from "store/slice/Workflow/WorkflowActions"
import { getWorkspace } from "store/slice/Workspace/WorkspaceActions"

/**
 * Get the appropriate file type for the initial node
 */
const getInitialNodeFileType = (_workspaceType?: number) => {
  return FILE_TYPE_SET.IMAGE
}

/**
 * Create initial state with workspace type
 * @param workspaceType - Workspace type to use directly (optional)
 * @param preserveFrom - Old state to preserve workspace type from (optional)
 * @param includeInitialNode - Whether to include initial node (default: true)
 */
const createInitialState = (
  workspaceType?: number,
  preserveFrom?: InputNode,
  includeInitialNode = true,
): InputNode => {
  const effectiveWorkspaceType =
    workspaceType !== undefined
      ? workspaceType
      : preserveFrom?.[WORKSPACE_TYPE_KEY]

  const newState: InputNode = {}

  // Add initial node if requested
  if (includeInitialNode) {
    newState[INITIAL_IMAGE_ELEMENT_ID] = {
      fileType: getInitialNodeFileType(effectiveWorkspaceType),
      param: {},
    } as InputNodeType
  }

  // Preserve workspace type in new state
  if (effectiveWorkspaceType !== undefined) {
    newState[WORKSPACE_TYPE_KEY] = effectiveWorkspaceType
  }

  return newState
}

const isPath = (p: unknown): p is string => typeof p === "string" && p !== ""

const normalizeFilePath = (
  path: unknown,
  filePathType: "single" | "array",
): string | string[] | undefined => {
  const paths = (Array.isArray(path) ? path : [path]).filter(isPath)
  if (filePathType === "array") {
    return paths.length > 0 ? paths : undefined
  }
  return paths[0]
}

const initialState: InputNode = createInitialState()

export const inputNodeSlice = createSlice({
  name: INPUT_NODE_SLICE_NAME,
  initialState,
  reducers: {
    deleteInputNode(state, action: PayloadAction<string>) {
      delete state[action.payload]
    },
    setCsvInputNodeParam(
      state,
      action: PayloadAction<{
        nodeId: string
        param: CsvInputParamType
      }>,
    ) {
      const { nodeId, param } = action.payload
      const inputNode = state[nodeId]
      if (isCsvInputNode(inputNode)) {
        inputNode.param = param
      }
    },
    setInputNodeMatlabPath(
      state,
      action: PayloadAction<{
        nodeId: string
        path: string
      }>,
    ) {
      const { nodeId, path } = action.payload
      const item = state[nodeId]
      if (isMatlabInputNode(item)) {
        item.matPath = path
      }
    },
    setInputNodeHDF5Path(
      state,
      action: PayloadAction<{
        nodeId: string
        path: string
      }>,
    ) {
      const { nodeId, path } = action.payload
      const item = state[nodeId]
      if (isHDF5InputNode(item)) {
        item.hdf5Path = path
      }
    },
  },
  extraReducers: (builder) =>
    builder
      .addCase(getWorkspace.fulfilled, (state) => {
        // Save workspace type for later use
        const workspaceType = WORKSPACE_TYPE.DEFAULT // Currently a fixed value
        state[WORKSPACE_TYPE_KEY] = workspaceType

        // Update the initial node based on workspace type
        const initialNode = state[INITIAL_IMAGE_ELEMENT_ID]
        if (initialNode) {
          const newFileType = getInitialNodeFileType(workspaceType)
          initialNode.fileType = newFileType as typeof initialNode.fileType
        }
      })
      .addCase(setInputNodeFilePath, (state, action) => {
        const { nodeId, filePath } = action.payload
        const targetNode = state[nodeId]
        targetNode.selectedFilePath = filePath
        if (isHDF5InputNode(targetNode)) {
          targetNode.hdf5Path = undefined
        }
        if (isMatlabInputNode(targetNode)) {
          targetNode.matPath = undefined
        }
      })
      .addCase(addInputNode, (state, action) => {
        const { node, fileType } = action.payload
        if (node.data?.type === NODE_TYPE_SET.INPUT) {
          try {
            const inputNode = FileNodeFactory.createInputNode(fileType)
            state[node.id] = inputNode
          } catch (error) {
            // eslint-disable-next-line no-console
            console.warn(`Unsupported file type: ${fileType}`, error)
          }
        }
      })
      .addCase(clearFlowElements, (state, action) =>
        createInitialState(action.payload?.workspaceType, state),
      )
      .addCase(deleteFlowNodes, (state, action) => {
        action.payload.forEach((node) => {
          if (node.data?.type === NODE_TYPE_SET.INPUT) {
            delete state[node.id]
          }
        })
      })
      .addCase(deleteFlowNodeById, (state, action) => {
        if (Object.keys(state).includes(action.payload)) {
          delete state[action.payload]
        }
      })
      .addCase(uploadFile.fulfilled, (state, action) => {
        const { nodeId } = action.meta.arg
        if (nodeId != null) {
          const { resultPath } = action.payload
          const target = state[nodeId]
          if (target) {
            const filePathType = FileNodeFactory.getFilePathType(
              target.fileType,
            )
            if (filePathType === "array") {
              target.selectedFilePath = [resultPath]
            } else {
              target.selectedFilePath = resultPath
            }
          }
        }
      })
      .addCase(fetchWorkflow.rejected, (state) => {
        // Clear previous workspace nodes and reset to initial state
        return createInitialState(undefined, state)
      })
      .addMatcher(
        isAnyOf(
          importWorkflowConfig.fulfilled,
          fetchWorkflow.fulfilled,
          reproduceWorkflow.fulfilled,
        ),
        (state, action) => {
          const newState = createInitialState(undefined, state, false)

          Object.values(action.payload.nodeDict)
            .filter(isInputNodePostData)
            .forEach((node) => {
              if (node.data?.fileType != null) {
                try {
                  const baseNode = FileNodeFactory.createInputNode(
                    node.data.fileType,
                  )
                  const filePathType = FileNodeFactory.getFilePathType(
                    node.data.fileType,
                  )
                  const specialPath = FileNodeFactory.getSpecialPathConfig(
                    node.data.fileType,
                  )

                  // Use specific param for CSV nodes
                  const param =
                    node.data.fileType === FILE_TYPE_SET.CSV
                      ? (node.data.param as CsvInputParamType)
                      : baseNode.param

                  const nodeState: InputNodeType = {
                    ...baseNode,
                    param,
                    selectedFilePath: normalizeFilePath(
                      node.data.path,
                      filePathType,
                    ),
                  } as InputNodeType

                  // Add special path properties if configured
                  if (specialPath) {
                    if (
                      specialPath.type === "hdf5Path" &&
                      "hdf5Path" in node.data &&
                      isHDF5InputNode(nodeState)
                    ) {
                      ;(nodeState as HDF5InputNode).hdf5Path = isPath(
                        node.data.hdf5Path,
                      )
                        ? node.data.hdf5Path
                        : undefined
                    } else if (
                      specialPath.type === "matPath" &&
                      "matPath" in node.data &&
                      isMatlabInputNode(nodeState)
                    ) {
                      ;(nodeState as MatlabInputNode).matPath = isPath(
                        node.data.matPath,
                      )
                        ? node.data.matPath
                        : undefined
                    }
                  }

                  newState[node.id] = nodeState
                } catch (error) {
                  // eslint-disable-next-line no-console
                  console.warn(
                    `Unsupported file type: ${node.data.fileType}`,
                    error,
                  )
                }
              }
            })
          return newState
        },
      ),
})

export const {
  setCsvInputNodeParam,
  setInputNodeMatlabPath,
  setInputNodeHDF5Path,
} = inputNodeSlice.actions

export default inputNodeSlice.reducer
