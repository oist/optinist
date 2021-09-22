import { createSlice, PayloadAction } from '@reduxjs/toolkit'
import { Elements } from 'react-flow-renderer'
import { initialElements } from 'const/flowchart'

export interface Param {
  name: string
  value: number
}
export interface Algorithm {
  name: string
  params: Array<Param>
}

export type Param_ = {
  [name: string]: number
}

export type Algorithm_ = {
  [name: string]: Param_
}
export interface Element {
  flowElements: Elements
  currentElement: string
  algoParams: Algorithm_
}

const initialAlgoParams_sample = {
  caiman_mc: {
    alpha_caiman: 30,
    beta_caiman: 30,
  },
  caiman_cnmf: {
    alpha_suite2p: 30,
    beta_suite2p: 30,
  },
}

const initialState: Element = {
  flowElements: initialElements,
  currentElement: 'caiman_mc',
  algoParams: initialAlgoParams_sample,
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
    setCurrentElement: (state, action: PayloadAction<string>) => {
      state.currentElement = action.payload
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
})

export const {
  setFlowElements,
  setElementPath,
  setCurrentElement,
  updateParam,
} = elementSlice.actions

export default elementSlice.reducer
