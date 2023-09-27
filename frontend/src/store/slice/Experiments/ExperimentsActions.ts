import { createAsyncThunk } from '@reduxjs/toolkit'
import {
  ExperimentsDTO,
  FetchExperimentDTO,
  getExperimentsApi,
  deleteExperimentByUidApi,
  deleteExperimentByListApi,
  fetchExperimentApi,
} from 'api/experiments/Experiments'
import { EXPERIMENTS_SLICE_NAME } from './ExperimentsType'
import { selectCurrentWorkspaceId } from '../Workspace/WorkspaceSelector'
import { ThunkApiConfig } from 'store/store'

export const getExperiments = createAsyncThunk<
  ExperimentsDTO,
  undefined,
  ThunkApiConfig
>(`${EXPERIMENTS_SLICE_NAME}/getExperiments`, async (_, thunkAPI) => {
  const workspaceId = selectCurrentWorkspaceId(thunkAPI.getState())
  if (workspaceId) {
    try {
      const response = await getExperimentsApi(workspaceId)
      return response
    } catch (e) {
      return thunkAPI.rejectWithValue(e)
    }
  } else {
    return thunkAPI.rejectWithValue('workspace id does not exist.')
  }
})

export const deleteExperimentByUid = createAsyncThunk<
  boolean,
  string,
  ThunkApiConfig
>(`${EXPERIMENTS_SLICE_NAME}/deleteExperimentByUid`, async (uid, thunkAPI) => {
  const workspaceId = selectCurrentWorkspaceId(thunkAPI.getState())
  if (workspaceId) {
    try {
      const response = await deleteExperimentByUidApi(workspaceId, uid)
      return response
    } catch (e) {
      return thunkAPI.rejectWithValue(e)
    }
  } else {
    return thunkAPI.rejectWithValue('workspace id does not exist.')
  }
})

export const deleteExperimentByList = createAsyncThunk<
  boolean,
  string[],
  ThunkApiConfig
>(`${EXPERIMENTS_SLICE_NAME}/deleteExperimentByList`, async (uid, thunkAPI) => {
  const workspaceId = selectCurrentWorkspaceId(thunkAPI.getState())
  if (workspaceId) {
    try {
      const response = await deleteExperimentByListApi(workspaceId, uid)
      return response
    } catch (e) {
      return thunkAPI.rejectWithValue(e)
    }
  } else {
    return thunkAPI.rejectWithValue('workspace id does not exist.')
  }
})

export const fetchExperiment = createAsyncThunk<FetchExperimentDTO, number>(
  `${EXPERIMENTS_SLICE_NAME}/fetchExperiment`,
  async (workspaceId, thunkAPI) => {
    try {
      const response = await fetchExperimentApi(workspaceId)
      return response
    } catch (e) {
      return thunkAPI.rejectWithValue(e)
    }
  },
)
