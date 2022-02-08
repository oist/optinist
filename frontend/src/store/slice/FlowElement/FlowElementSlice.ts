import { createSlice, PayloadAction } from '@reduxjs/toolkit'
import {
  Elements,
  removeElements,
  Node,
  Position,
  isNode,
  FlowTransform,
} from 'react-flow-renderer'
import {
  FLOW_ELEMENT_SLICE_NAME,
  FlowElement,
  NODE_TYPE_SET,
  NodeData,
  ElementCoord,
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

const initialFlowPosition: FlowTransform = {
  x: 0,
  y: 0,
  zoom: 0.8,
}

const initialElementCoord: ElementCoord = {
  x: 350,
  y: 60,
}

const initialState: FlowElement = {
  flowElements: initialElements,
  flowPosition: initialFlowPosition,
  elementCoord: initialElementCoord,
}

export const flowElementSlice = createSlice({
  name: FLOW_ELEMENT_SLICE_NAME,
  initialState,
  reducers: {
    setFlowPosition: (state, action: PayloadAction<FlowTransform>) => {
      state.flowPosition = action.payload
    },
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
    edifFlowElementsLabelById: (
      state,
      action: PayloadAction<{
        nodeId: string
        fileName: string
      }>,
    ) => {
      let { nodeId, fileName } = action.payload
      const elementIdx = state.flowElements.findIndex(
        (ele) => ele.id === nodeId,
      )
      if (state.flowElements[elementIdx].data?.label) {
        state.flowElements[elementIdx].data!.label = fileName
      }
    },
    editFlowElementPositionById: (
      state,
      action: PayloadAction<{
        nodeId: string
        coord: {
          x: number
          y: number
        }
      }>,
    ) => {
      let { nodeId, coord } = action.payload
      const elementIdx = state.flowElements.findIndex(
        (ele) => ele.id === nodeId,
      )
      const targetItem = state.flowElements[elementIdx]
      if (isNode(targetItem)) {
        targetItem.position = coord
      }
    },
    setElementCoord: (state) => {
      const x = state.elementCoord.x
      const y = state.elementCoord.y
      if (x > 800 || y > 200) {
        state.elementCoord.x = 300
        state.elementCoord.y = 100
      } else {
        state.elementCoord.x += 250
      }
    },
  },
})

export const {
  setFlowPosition,
  setFlowElements,
  addFlowElementNode,
  deleteFlowElements,
  deleteFlowElementsById,
  edifFlowElementsLabelById,
  editFlowElementPositionById,
  setElementCoord,
} = flowElementSlice.actions

export default flowElementSlice.reducer
