import { createSlice, PayloadAction } from '@reduxjs/toolkit'
import { Elements, Node } from 'react-flow-renderer'
import { initialElements, INITIAL_ALGO_ELEMENT_ID } from 'const/flowchart'
import { getAlgoParams } from './ElementAction'
import { Element, NodeData } from './ElementType'

const initialState: Element = {
  flowElements: initialElements,
  currentElementId: INITIAL_ALGO_ELEMENT_ID,
  algoParams: {},
}

export const elementSlice = createSlice({
  name: 'element',
  initialState,
  reducers: {
    setFlowElements: (state, action: PayloadAction<Elements>) => {
      state.flowElements = action.payload
    },
    addFlowElement: (state, action: PayloadAction<Node<NodeData>>) => {
      state.flowElements.push(action.payload)
    },
    setElementPath: (
      state,
      action: PayloadAction<{ id: string; path: string }>,
    ) => {
      var idx = state.flowElements.findIndex((e) => e.id === action.payload.id)
      const node = state.flowElements[idx]
      if (node && node.data) {
        node.data.path = action.payload.path
      }
    },
    setCurrentElement: (state, action: PayloadAction<string>) => {
      state.currentElementId = action.payload
    },
    updateParam: (
      state,
      action: PayloadAction<{ paramKey: string; newValue: unknown }>,
    ) => {
      const { paramKey, newValue } = action.payload
      state.algoParams[state.currentElementId].param[paramKey] = newValue
    },
  },
  extraReducers: (builder) => {
    builder.addCase(getAlgoParams.fulfilled, (state, action) => {
      const { id, algoName } = action.meta.arg
      state.algoParams[id] = {
        name: algoName,
        param: action.payload,
      }
    })
  },
})

export const {
  setFlowElements,
  setElementPath,
  setCurrentElement,
  updateParam,
  addFlowElement,
} = elementSlice.actions

export default elementSlice.reducer
