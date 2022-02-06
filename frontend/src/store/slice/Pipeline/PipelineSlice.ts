import { createSlice, PayloadAction } from '@reduxjs/toolkit'
import { pollRunResult, run } from './PilepineActions'

import { Pipeline, PIPELINE_SLICE_NAME } from './PipelineType'

const initialState: Pipeline = {
  pipelines: {},
}

export const pipelineSlice = createSlice({
  name: PIPELINE_SLICE_NAME,
  initialState,
  reducers: {},
  extraReducers: (builder) => {
    builder
      .addCase(run.pending, (state, action) => {})
      .addCase(run.fulfilled, (state, action) => {})
      .addCase(run.rejected, (state, action) => {})
      .addCase(pollRunResult.pending, (state, action) => {})
      .addCase(pollRunResult.fulfilled, (state, action) => {})
      .addCase(pollRunResult.rejected, (state, action) => {})
  },
})

export const {} = pipelineSlice.actions

export default pipelineSlice.reducer
