import { createAsyncThunk } from '@reduxjs/toolkit'
import {
  ExperimentDTO,
  ExperimentsDTO,
  getExperimentsApi,
  deleteExperimentByUidApi,
  importExperimentByUidApi,
  deleteExperimentByListApi,
  fetchExperimentApi,
} from 'api/experiments/Experiments'
import { RunPostData } from 'api/run/Run'
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

export const fetchExperiment = createAsyncThunk<ExperimentDTO, undefined>(
  `${EXPERIMENTS_SLICE_NAME}/fetchExperiment`,
  async (_, thunkAPI) => {
    try {
      const response = await fetchExperimentApi()
      return (
        response ??
        thunkAPI.rejectWithValue({ message: 'No experiments found' })
      )
    } catch (e) {
      return thunkAPI.rejectWithValue(e)
    }
  },
)

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

export const deleteExperimentByList = createAsyncThunk<boolean, string[]>(
  `${EXPERIMENTS_SLICE_NAME}/deleteExperimentByList`,
  async (uid, thunkAPI) => {
    try {
      const response = await deleteExperimentByListApi(uid)
      return response
    } catch (e) {
      return thunkAPI.rejectWithValue(e)
    }
  },
)

export const importExperimentByUid = createAsyncThunk<RunPostData, string>(
  `${EXPERIMENTS_SLICE_NAME}/importExperimentByUid`,
  async (uid, thunkAPI) => {
    try {
      const response = await importExperimentByUidApi(uid)
      return response
    } catch (e) {
      return thunkAPI.rejectWithValue(e)
    }
  },
)
