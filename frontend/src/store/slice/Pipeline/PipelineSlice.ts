import { createSlice, isAnyOf, PayloadAction } from '@reduxjs/toolkit'
import { importExperimentByUid } from '../Experiments/ExperimentsActions'
import { pollRunResult, run, runByCurrentUid } from './PipelineActions'
import {
  Pipeline,
  PIPELINE_SLICE_NAME,
  RUN_BTN_OPTIONS,
  RUN_BTN_TYPE,
  RUN_STATUS,
} from './PipelineType'

import {
  getInitialRunResult,
  convertToRunResult,
  isNodeResultPending,
} from './PipelineUtils'

const initialState: Pipeline = {
  run: {
    status: RUN_STATUS.START_UNINITIALIZED,
  },
  runBtn: RUN_BTN_OPTIONS.RUN_NEW,
}

export const pipelineSlice = createSlice({
  name: PIPELINE_SLICE_NAME,
  initialState,
  reducers: {
    cancelPipeline(state) {
      state.run.status = RUN_STATUS.CANCELED
    },
    setRunBtnOption: (
      state,
      action: PayloadAction<{
        runBtnOption: RUN_BTN_TYPE
      }>,
    ) => {
      state.runBtn = action.payload.runBtnOption
    },
    clearCurrentPipeline: () => initialState,
  },
  extraReducers: (builder) => {
    builder
      .addCase(pollRunResult.fulfilled, (state, action) => {
        if (state.run.status === RUN_STATUS.START_SUCCESS) {
          state.run.runResult = {
            ...state.run.runResult, // pendingのNodeResultはそのままでsuccessもしくはerrorのみ上書き
            ...convertToRunResult(action.payload),
          }
          const runResultPendingList = Object.values(
            state.run.runResult,
          ).filter(isNodeResultPending)
          if (runResultPendingList.length === 0) {
            // 終了
            state.run.status = RUN_STATUS.FINISHED
          }
        }
      })
      .addCase(pollRunResult.rejected, (state, action) => {
        state.run.status = RUN_STATUS.ABORTED
      })
      .addCase(importExperimentByUid.fulfilled, (state, action) => {
        state.currentPipeline = {
          uid: action.meta.arg.uid,
        }
        state.runBtn = RUN_BTN_OPTIONS.RUN_ALREADY
        state.run = {
          status: RUN_STATUS.START_UNINITIALIZED,
        }
      })
      .addMatcher(
        isAnyOf(run.pending, runByCurrentUid.pending),
        (state, action) => {
          state.run = {
            status: RUN_STATUS.START_PENDING,
          }
        },
      )
      .addMatcher(
        isAnyOf(run.fulfilled, runByCurrentUid.fulfilled),
        (state, action) => {
          const runPostData = action.meta.arg.runPostData
          const uid = action.payload
          state.run = {
            uid,
            status: RUN_STATUS.START_SUCCESS,
            runResult: getInitialRunResult({ name: '', ...runPostData }),
            runPostData: { name: '', ...runPostData },
          }
          state.currentPipeline = {
            uid: action.payload,
          }
        },
      )
      .addMatcher(
        isAnyOf(run.rejected, runByCurrentUid.rejected),
        (state, action) => {
          state.run = {
            status: RUN_STATUS.START_ERROR,
          }
        },
      )
  },
})

export const { cancelPipeline, setRunBtnOption, clearCurrentPipeline } =
  pipelineSlice.actions

export default pipelineSlice.reducer
