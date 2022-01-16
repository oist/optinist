import { createSlice, PayloadAction } from '@reduxjs/toolkit'
import {
  Elements,
  removeElements,
  Node,
  Position,
  isEdge,
} from 'react-flow-renderer'
import {
  FLOW_ELEMENT_SLICE_NAME,
  FlowElement,
  NODE_TYPE_SET,
  NodeData,
} from './FlowElementType'
import {
  INITIAL_ALGO_STYLE,
  INITIAL_DATA_STYLE,
  INITIAL_IMAGE_ELEMENT_ID,
  INITIAL_IMAGE_ELEMENT_NAME,
} from 'const/flowchart'
import { FILE_TYPE } from '../InputNode/InputNodeType'

const initialElements: Elements<NodeData> = [
  {
    id: INITIAL_IMAGE_ELEMENT_ID,
    type: 'ImageFileNode',
    data: {
      type: NODE_TYPE_SET.INPUT,
      label: INITIAL_IMAGE_ELEMENT_NAME,
    },
    style: INITIAL_DATA_STYLE,
    position: { x: 50, y: 60 },
  },
]

const initialState: FlowElement = {
  flowElements: initialElements,
}

export const flowElementSlice = createSlice({
  name: FLOW_ELEMENT_SLICE_NAME,
  initialState,
  reducers: {
    setFlowElements: (state, action: PayloadAction<Elements>) => {
      state.flowElements = action.payload
    },
    deleteFlowElements: (state, action: PayloadAction<Elements>) => {
      state.flowElements = removeElements(action.payload, state.flowElements)
    },
    deleteFlowElementsById: (state, action: PayloadAction<string>) => {
      const element = state.flowElements.find(
        (edge) => edge.id === action.payload,
      )
      if (element !== undefined) {
        state.flowElements = removeElements([element], state.flowElements)
      }
    },
    addFlowElementNode: (
      state,
      action: PayloadAction<{
        node: Node<NodeData>
        inputNodeInfo?: { fileType: FILE_TYPE }
        algoNodeInfo?: { functionPath: string; name: string }
      }>,
    ) => {
      let { node } = action.payload
      if (node.data?.type === NODE_TYPE_SET.INPUT) {
        node = {
          ...node,
          style: {
            ...node.style,
            ...INITIAL_DATA_STYLE,
          },
          targetPosition: Position.Left,
          sourcePosition: Position.Right,
        }
      } else if (node.data?.type === NODE_TYPE_SET.ALGORITHM) {
        node = {
          ...node,
          style: {
            ...node.style,
            ...INITIAL_ALGO_STYLE,
          },
          targetPosition: Position.Left,
          sourcePosition: Position.Right,
        }
      }
      state.flowElements.push(node)
    },
  },
})

export const {
  setFlowElements,
  addFlowElementNode,
  deleteFlowElements,
  deleteFlowElementsById,
} = flowElementSlice.actions

export default flowElementSlice.reducer
