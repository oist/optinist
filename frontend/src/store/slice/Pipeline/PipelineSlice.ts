import { createSlice, PayloadAction } from '@reduxjs/toolkit'
import { pollRunResult, run } from './PipelineActions'
import { Pipeline, PIPELINE_SLICE_NAME, RUN_STATUS } from './PipelineType'

import {
  getInitialRunResult,
  convertToRunResult,
  isStartedPipeline,
  isNodeResultPending,
} from './PipelineUtils'

const initialState: Pipeline = {
  pipelines: {},
  uidHistory: [],
}

export const pipelineSlice = createSlice({
  name: PIPELINE_SLICE_NAME,
  initialState,
  reducers: {
    cancelPipeline(state, action: PayloadAction<string>) {
      const target = Object.values(state.pipelines).find(
        (pipeline) =>
          pipeline.status === RUN_STATUS.START_SUCCESS &&
          pipeline.uid === action.payload,
      )
      if (target != null) {
        target.status = RUN_STATUS.CANCELED
        state.uidHistory.pop()
      }
    },
  },
  extraReducers: (builder) => {
    builder
      .addCase(run.pending, (state, action) => {
        state.pipelines[action.meta.requestId] = {
          status: RUN_STATUS.START_PENDING,
        }
      })
      .addCase(run.fulfilled, (state, action) => {
        const runPostData = action.meta.arg.runPostData
        const uid = action.payload
        state.pipelines[action.meta.requestId] = {
          uid,
          status: RUN_STATUS.START_SUCCESS,
          runResult: getInitialRunResult(runPostData),
          runPostData,
        }
        state.uidHistory.push(uid)
      })
      .addCase(run.rejected, (state, action) => {
        state.pipelines[action.meta.requestId] = {
          status: RUN_STATUS.START_ERROR,
        }
      })
      .addCase(pollRunResult.fulfilled, (state, action) => {
        const target = Object.values(state.pipelines).find(
          (pipeline) =>
            pipeline.status === RUN_STATUS.START_SUCCESS &&
            pipeline.uid === action.meta.arg.uid,
        )
        if (target != null && target.status === RUN_STATUS.START_SUCCESS) {
          target.runResult = {
            ...target.runResult, // pendingのNodeResultはそのままでsuccessもしくはerrorのみ上書き
            ...convertToRunResult(action.payload),
          }
          const runResultPendingList = Object.values(target.runResult).filter(
            isNodeResultPending,
          )
          if (runResultPendingList.length === 0) {
            // 終了
            target.status = RUN_STATUS.FINISHED
          }
        }
      })
      .addCase(pollRunResult.rejected, (state, action) => {
        const target = Object.values(state.pipelines).find(
          (pipeline) =>
            isStartedPipeline(pipeline) && pipeline.uid === action.meta.arg.uid,
        )
        if (target != null) {
          target.status = RUN_STATUS.ABORTED
        }
      })
  },
})

export const { cancelPipeline } = pipelineSlice.actions

export default pipelineSlice.reducer
