import { createSlice, PayloadAction } from '@reduxjs/toolkit'
import { Elements, Node, Position } from 'react-flow-renderer'
import {
  initialElements,
  INITIAL_ALGO_STYLE,
  INITIAL_DATA_STYLE,
} from 'const/flowchart'
import { NodeData, NODE_DATA_TYPE_SET } from 'const/NodeData'
import {
  isAlgoNodeData,
  isCsvNodeData,
  isImageNodeData,
} from 'utils/ElementUtils'
import { Element, ELEMENT_SLICE_NAME } from './ElementType'
import { uploadCsvFile, uploadImageFile } from '../FileData/FileDataAction'
import { selectCsvFile, selectImageFile } from '../FileData/FileData'

const initialState: Element = {
  flowElements: initialElements,
}

export const elementSlice = createSlice({
  name: ELEMENT_SLICE_NAME,
  initialState,
  reducers: {
    setFlowElements: (state, action: PayloadAction<Elements>) => {
      state.flowElements = action.payload
    },
    addFlowElement: (state, action: PayloadAction<Node<NodeData>>) => {
      let node = action.payload
      if (isAlgoNodeData(node)) {
        node = {
          ...node,
          style: {
            ...node.style,
            ...INITIAL_ALGO_STYLE,
          },
          targetPosition: Position.Left,
          sourcePosition: Position.Right,
        }
      } else if (isImageNodeData(node) || isCsvNodeData(node)) {
        node = {
          ...node,
          style: {
            ...node.style,
            ...INITIAL_DATA_STYLE,
          },
          targetPosition: Position.Left,
          sourcePosition: Position.Right,
        }
      }
      state.flowElements.push(node)
    },
  },
  extraReducers: (builder) => {
    builder
      .addCase(uploadImageFile.fulfilled, (state, action) => {
        const { nodeId } = action.meta.arg
        const { path } = action.payload
        setFilePathFn(state, nodeId, path)
      })
      .addCase(uploadImageFile.pending, (state, action) => {
        const { nodeId } = action.meta.arg
        setFilePathFn(state, nodeId, undefined)
      })
      .addCase(selectImageFile, (state, action) => {
        const { nodeId, path } = action.payload
        setFilePathFn(state, nodeId, path)
      })
      .addCase(uploadCsvFile.fulfilled, (state, action) => {
        const { nodeId } = action.meta.arg
        const { path } = action.payload
        setFilePathFn(state, nodeId, path)
      })
      .addCase(uploadCsvFile.pending, (state, action) => {
        const { nodeId } = action.meta.arg
        setFilePathFn(state, nodeId, undefined)
      })
      .addCase(selectCsvFile, (state, action) => {
        const { nodeId, path } = action.payload
        setFilePathFn(state, nodeId, path)
      })
  },
})

function setFilePathFn(
  state: Element,
  nodeId: string,
  path: string | undefined,
) {
  const idx = state.flowElements.findIndex((e) => e.id === nodeId)
  const node = state.flowElements[idx]
  if (isImageNodeData(node) && node.data) {
    node.data = {
      ...node.data,
      type: NODE_DATA_TYPE_SET.IMAGE,
      path,
    }
  } else if (isCsvNodeData(node) && node.data) {
    node.data = {
      ...node.data,
      type: NODE_DATA_TYPE_SET.CSV,
      path,
    }
  }
}

export const { setFlowElements, addFlowElement } = elementSlice.actions

export default elementSlice.reducer
