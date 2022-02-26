import { createSlice } from '@reduxjs/toolkit'
import { EXPERIMENTS_SLICE_NAME, Experiments } from './ExperimentsType'
import {
  getExperiments,
  //   getExperimentByUid,
  deleteExperimentByUid,
} from './ExperimentsActions'
import {
  convertToExperimentListType,
  //   convertToExperimentType,
} from './ExperimentsUtils'
import { run } from '../Pipeline/PipelineActions'

const initialState: Experiments = {
  status: 'uninitialized',
}

export const experimentsSlice = createSlice({
  name: EXPERIMENTS_SLICE_NAME,
  initialState: initialState as Experiments,
  reducers: {},
  extraReducers: (builder) => {
    builder
      .addCase(run.fulfilled, () => {
        return {
          status: 'uninitialized',
        }
      })
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
  },
})

export default experimentsSlice.reducer
