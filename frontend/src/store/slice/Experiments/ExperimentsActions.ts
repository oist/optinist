import { createAsyncThunk } from '@reduxjs/toolkit'
import {
  ExperimentsDTO,
  //   ExperimentDTO,
  getExperiments as getExperimentsApi,
  //   getExperimentByUid as getExperimentByUidApi,
  deleteExperimentByUid as deleteExperimentByUidApi,
} from 'api/experiments/Experiments'
import { EXPERIMENTS_SLICE_NAME } from './ExperimentsType'

export const getExperiments = createAsyncThunk<ExperimentsDTO, undefined>(
  `${EXPERIMENTS_SLICE_NAME}/getExperiments`,
  async (_, thunkAPI) => {
    try {
      const response = await getExperimentsApi()
      return response
    } catch (e) {
      return thunkAPI.rejectWithValue(e)
    }
  },
)

// export const getExperimentByUid = createAsyncThunk<{}, string>(
//   `${EXPERIMENTS_SLICE_NAME}/getExperimentByUid`,
//   async (uid, thunkAPI) => {
//     try {
//       const response = await getExperimentByUidApi(uid)
//       return response
//     } catch (e) {
//       return thunkAPI.rejectWithValue(e)
//     }
//   },
// )

export const deleteExperimentByUid = createAsyncThunk<boolean, string>(
  `${EXPERIMENTS_SLICE_NAME}/deleteExperimentByUid`,
  async (uid, thunkAPI) => {
    try {
      const response = await deleteExperimentByUidApi(uid)
      return response
    } catch (e) {
      return thunkAPI.rejectWithValue(e)
    }
  },
)
