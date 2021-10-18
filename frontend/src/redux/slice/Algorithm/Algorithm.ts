import { createSlice, PayloadAction } from '@reduxjs/toolkit'

import { INITIAL_ALGO_ELEMENT_ID } from 'const/flowchart'
import { NODE_DATA_TYPE_SET } from 'const/NodeData'
import { clickNode, runPipeline } from '../Element/ElementAction'
import { getAlgoOutputData, getAlgoParams } from './AlgorithmAction'

import { ALGORITHM_SLICE_NAME, Algorithm } from './AlgorithmType'

const initialState: Algorithm = {
  currentAlgoId: INITIAL_ALGO_ELEMENT_ID,
  algoMap: {},
}

export const algorithmSlice = createSlice({
  name: ALGORITHM_SLICE_NAME,
  initialState,
  reducers: {
    updateParam: (
      state,
      action: PayloadAction<{ paramKey: string; newValue: unknown }>,
    ) => {
      const { paramKey, newValue } = action.payload
      state.algoMap[state.currentAlgoId].param[paramKey] = newValue
    },
  },
  extraReducers: (builder) => {
    builder
      .addCase(clickNode, (state, action) => {
        if (action.payload.type === NODE_DATA_TYPE_SET.ALGO) {
          state.currentAlgoId = action.payload.id
        }
      })
      .addCase(getAlgoParams.fulfilled, (state, action) => {
        const { id, algoName } = action.meta.arg
        const algo = state.algoMap[id]
        state.algoMap[id] = {
          ...algo,
          name: algoName,
          param: action.payload,
        }
      })
      // .addCase(getAlgoOutputData.fulfilled, (state, action) => {
      //   console.log(action.payload)
      //   state.algoMap[action.meta.arg.id].output = { data: action.payload }
      // })
      .addCase(runPipeline.fulfilled, (state, action) => {
        if (action.payload.message === 'success') {
          console.log(action.payload.outputPaths)
          // state.algoMap[state.currentAlgoId].output = {
          //   imageDir: action.payload.outputPaths,
          // }
        }
      })
  },
})

export const { updateParam } = algorithmSlice.actions

export default algorithmSlice.reducer
