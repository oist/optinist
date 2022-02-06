import { createAsyncThunk } from '@reduxjs/toolkit'

import { ThunkApiConfig } from 'store/store'
import { selectElementListForRun } from '../FlowElement/FlowElementSelectors'
import { selectNwbParams } from '../NWB/NWBSelectors'
import { selectSnakemakeParams } from '../Snakemake/SnakemakeSelectors'
import { PIPELINE_SLICE_NAME } from './PipelineType'
import {
  run as runRequest,
  runResult,
  RunResultDTO,
  RunPostData,
} from 'api/Run/Run'

export const run = createAsyncThunk<string, string | undefined, ThunkApiConfig>(
  `${PIPELINE_SLICE_NAME}/run`,
  async (uid, thunkAPI) => {
    const nwbParam = selectNwbParams(thunkAPI.getState())
    const snakemakeParam = selectSnakemakeParams(thunkAPI.getState())
    const elementListForRun = selectElementListForRun(thunkAPI.getState())
    const runPostData: RunPostData = {
      nwbParam,
      snakemakeParam,
      ...elementListForRun,
    }
    try {
      const responseData = await runRequest({ runData: runPostData, uid })
      return responseData
    } catch (e) {
      return thunkAPI.rejectWithValue(e)
    }
  },
)

export const pollRunResult = createAsyncThunk<
  RunResultDTO,
  {
    //   todo
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
