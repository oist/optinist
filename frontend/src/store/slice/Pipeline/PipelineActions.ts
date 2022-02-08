import { createAsyncThunk } from '@reduxjs/toolkit'

import { ThunkApiConfig } from 'store/store'
import { PIPELINE_SLICE_NAME } from './PipelineType'
import {
  run as runRequest,
  runResult,
  RunResultDTO,
  RunPostData,
} from 'api/Run/Run'

export const run = createAsyncThunk<
  string,
  { uid: string | undefined; runPostData: RunPostData },
  ThunkApiConfig
>(`${PIPELINE_SLICE_NAME}/run`, async ({ uid, runPostData }, thunkAPI) => {
  try {
    const responseData = await runRequest({ runData: runPostData, uid })
    return responseData
  } catch (e) {
    return thunkAPI.rejectWithValue(e)
  }
})

export const pollRunResult = createAsyncThunk<
  RunResultDTO,
  {
    uid: string
  },
  ThunkApiConfig
>(`${PIPELINE_SLICE_NAME}/pollRunResult`, async ({ uid }, thunkAPI) => {
  try {
    const responseData = await runResult({ uid })
    return responseData
  } catch (e) {
    return thunkAPI.rejectWithValue(e)
  }
})
