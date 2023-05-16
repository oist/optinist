import { createSlice, isAnyOf } from '@reduxjs/toolkit'
import { EXPERIMENTS_SLICE_NAME, Experiments } from './ExperimentsType'
import {
  getExperiments,
  deleteExperimentByUid,
  deleteExperimentByList,
  fetchExperiment,
} from './ExperimentsActions'
import {
  convertToExperimentListType,
  convertToExperimentType,
} from './ExperimentsUtils'
import {
  pollRunResult,
  run,
  runByCurrentUid,
} from '../Pipeline/PipelineActions'

export const initialState: Experiments = {
  status: 'uninitialized',
}

export const experimentsSlice = createSlice({
  name: EXPERIMENTS_SLICE_NAME,
  initialState: initialState as Experiments,
  reducers: {},
  extraReducers: (builder) => {
    builder
      .addCase(getExperiments.pending, () => {
        return {
          status: 'pending',
        }
      })
      .addCase(getExperiments.fulfilled, (state, action) => {
        const experimentList = convertToExperimentListType(action.payload)
        return {
          status: 'fulfilled',
          experimentList,
        }
      })
      .addCase(getExperiments.rejected, (state, action) => {
        return {
          status: 'error',
          message: action.error.message,
        }
      })
      .addCase(deleteExperimentByUid.fulfilled, (state, action) => {
        if (action.payload && state.status === 'fulfilled') {
          delete state.experimentList[action.meta.arg]
        }
      })
      .addCase(deleteExperimentByList.fulfilled, (state, action) => {
        if (action.payload && state.status === 'fulfilled') {
          action.meta.arg.map((v) => delete state.experimentList[v])
        }
      })
      .addCase(pollRunResult.fulfilled, (state, action) => {
        if (state.status === 'fulfilled') {
          const uid = action.meta.arg.uid
          const target = state.experimentList[uid]
          Object.entries(action.payload).forEach(([nodeId, value]) => {
            if (value.status === 'success') {
              target.functions[nodeId].status = 'success'
            } else if (value.status === 'error') {
              target.functions[nodeId].status = 'error'
            }
          })
        }
      })
      .addCase(fetchExperiment.fulfilled, (state, action) => {
        if (state.status === 'fulfilled') {
          state.experimentList[action.payload.unique_id] =
            convertToExperimentType(action.payload)
        }
      })
      .addMatcher(isAnyOf(run.fulfilled, runByCurrentUid.fulfilled), () => {
        return {
          status: 'uninitialized',
        }
      })
  },
})

export default experimentsSlice.reducer
