import { createAsyncThunk } from '@reduxjs/toolkit'

import { ThunkApiConfig } from 'store/store'
import { PIPELINE_SLICE_NAME } from './PipelineType'
import {
  runApi,
  runByUidApi,
  runResult,
  RunResultDTO,
  RunPostData,
} from 'api/run/Run'
import {
  selectPipelineLatestUid,
  selectRunResultPendingNodeIdList,
} from './PipelineSelectors'

export const run = createAsyncThunk<
  string,
  { runPostData: RunPostData },
  ThunkApiConfig
>(`${PIPELINE_SLICE_NAME}/run`, async ({ runPostData }, thunkAPI) => {
  try {
    const responseData = await runApi(runPostData)
    return responseData
  } catch (e) {
    return thunkAPI.rejectWithValue(e)
  }
})

export const runByCurrentUid = createAsyncThunk<
  string,
  { runPostData: Omit<RunPostData, 'name'> },
  ThunkApiConfig
>(`${PIPELINE_SLICE_NAME}/run`, async ({ runPostData }, thunkAPI) => {
  const currentUid = selectPipelineLatestUid(thunkAPI.getState())
  if (currentUid != null) {
    try {
      const responseData = await runByUidApi(currentUid, runPostData)
      return responseData
    } catch (e) {
      return thunkAPI.rejectWithValue(e)
    }
  } else {
    return thunkAPI.rejectWithValue('currentUid dose not exist.')
  }
})

export const pollRunResult = createAsyncThunk<
  RunResultDTO,
  {
    uid: string
  },
  ThunkApiConfig
>(`${PIPELINE_SLICE_NAME}/pollRunResult`, async ({ uid }, thunkAPI) => {
  const pendingNodeIdList = selectRunResultPendingNodeIdList(
    thunkAPI.getState(),
  )
  try {
    const responseData = await runResult({ uid, pendingNodeIdList })
    return responseData
  } catch (e) {
    return thunkAPI.rejectWithValue(e)
  }
})
