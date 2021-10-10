import { createSlice, PayloadAction } from '@reduxjs/toolkit'
import { INITIAL_OUTPUT_ELEMENT_ID } from 'const/flowchart'
import { clickNode, runPipeline } from '../Element/ElementAction'
import { getOutputData } from './OutputAction'
import { OUTPUT_SLICE_NAME, Output } from './OutputType'

const initialState: Output = {
  currentOutputId: INITIAL_OUTPUT_ELEMENT_ID,
  outputData: {},
}

export const outputSlice = createSlice({
  name: OUTPUT_SLICE_NAME,
  initialState,
  reducers: {},
  extraReducers: (builder) => {
    builder
      .addCase(getOutputData.fulfilled, (state, action) => {
        state.outputData[action.meta.arg.id] = action.payload
      })
      .addCase(clickNode, (state, action) => {
        if (action.payload.type === 'output') {
          state.currentOutputId = action.payload.id
        }
      })
      .addCase(runPipeline.fulfilled, (state, action) => {
        if (action.payload.message === 'success') {
          state.outputData[state.currentOutputId] = action.payload.data
        }
      })
  },
})

export const {} = outputSlice.actions

export default outputSlice.reducer
