import { createSlice, PayloadAction, current } from '@reduxjs/toolkit'
import { Elements } from 'react-flow-renderer'
import { initialElements } from 'const/flowchart'
import { getAlgoParams } from './ElementAction'
import { Element } from './ElementType'

const initialState: Element = {
  flowElements: initialElements,
  currentElement: 'caiman_mc',
  algoParams: {},
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
      var idx = state.flowElements.findIndex((e) => e.id === action.payload.id)
      state.flowElements[idx].data.path = action.payload.path
    },
    setCurrentElement: (state, action: PayloadAction<string>) => {
      state.currentElement = action.payload
      console.log(current(state))
    },
    updateParam: (
      state,
      action: PayloadAction<{ name: string; newValue: number }>,
    ) => {
      const { name, newValue } = action.payload
      state.algoParams[state.currentElement][name] = newValue
      console.log(state.algoParams[state.currentElement][name])
    },
  },
  extraReducers: (builder) => {
    // Add reducers for additional action types here, and handle loading state as needed
    builder.addCase(getAlgoParams.fulfilled, (state, action) => {
      // Add user to the state array
      state.algoParams[state.currentElement] = action.payload
      console.log(current(state))
    })
  },
})

export const {
  setFlowElements,
  setElementPath,
  setCurrentElement,
  updateParam,
} = elementSlice.actions

export default elementSlice.reducer
