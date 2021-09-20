import { createSlice, PayloadAction } from '@reduxjs/toolkit'
import { Elements } from 'react-flow-renderer'
import { initialElements } from 'const/flowchart'

export interface Element {
  flowElements: Elements
}

const initialState: Element = {
  flowElements: initialElements,
}

export const elementSlice = createSlice({
  name: 'element',
  initialState,
  // The `reducers` field lets us define reducers and generate associated actions
  reducers: {
    setFlowElements: (state, action: PayloadAction<Elements>) => {
      state.flowElements = action.payload
    },
    setElementPath: (
      state,
      action: PayloadAction<{ id: string; path: string }>,
    ) => {
      var idx = state.flowElements.findIndex((e) => e.id == action.payload.id)
      state.flowElements[idx].data.path = action.payload.path
    },
  },
})

export const { setFlowElements, setElementPath } = elementSlice.actions

export default elementSlice.reducer
