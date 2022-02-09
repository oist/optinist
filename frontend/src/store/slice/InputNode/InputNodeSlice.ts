import { createSlice, PayloadAction } from '@reduxjs/toolkit'
import { INITIAL_IMAGE_ELEMENT_ID } from 'const/flowchart'
import {
  addFlowElementNode,
  deleteFlowElements,
  deleteFlowElementsById,
} from '../FlowElement/FlowElementSlice'
import { NODE_TYPE_SET } from '../FlowElement/FlowElementType'
import { isNodeData } from '../FlowElement/FlowElementUtils'
import {
  FILE_TYPE_SET,
  InputNode,
  INPUT_NODE_SLICE_NAME,
} from './InputNodeType'

const initialState: InputNode = {
  [INITIAL_IMAGE_ELEMENT_ID]: {
    fileType: FILE_TYPE_SET.IMAGE,
  },
}

export const inputNodeSlice = createSlice({
  name: INPUT_NODE_SLICE_NAME,
  initialState,
  reducers: {
    deleteInputNode(state, action: PayloadAction<string>) {
      delete state[action.payload]
    },
    setInputImageNodeFile(
      state,
      action: PayloadAction<{
        nodeId: string
        filePath: string
      }>,
    ) {
      const { nodeId, filePath } = action.payload
      state[nodeId].selectedFilePath = filePath
    },
    setInputNodeFilePath(
      state,
      action: PayloadAction<{
        nodeId: string
        filePath: string
      }>,
    ) {
      const { nodeId, filePath } = action.payload
      state[nodeId].selectedFilePath = filePath
    },
  },
  extraReducers: (builder) =>
    builder
      .addCase(addFlowElementNode, (state, action) => {
        const { node, inputNodeInfo } = action.payload
        if (node.data?.type === NODE_TYPE_SET.INPUT && inputNodeInfo != null) {
          const fileType = inputNodeInfo.fileType
          switch (fileType) {
            case FILE_TYPE_SET.CSV:
              state[node.id] = {
                fileType,
              }
              break
            case FILE_TYPE_SET.IMAGE:
              state[node.id] = {
                fileType,
              }
              break
            case FILE_TYPE_SET.HDF5:
              state[node.id] = {
                fileType,
              }
              break
          }
        }
      })
      .addCase(deleteFlowElements, (state, action) => {
        action.payload
          .filter((node) => isNodeData(node))
          .forEach((node) => {
            if (node.data?.type === NODE_TYPE_SET.INPUT) {
              delete state[node.id]
            }
          })
      })
      .addCase(deleteFlowElementsById, (state, action) => {
        if (Object.keys(state).includes(action.payload)) {
          delete state[action.payload]
        }
      }),
})

export const { setInputNodeFilePath, setInputImageNodeFile } =
  inputNodeSlice.actions

export default inputNodeSlice.reducer
