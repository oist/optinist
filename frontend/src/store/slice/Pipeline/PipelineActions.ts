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
import { selectCurrentWorkspaceId } from '../Workspace/WorkspaceSelector'

export const run = createAsyncThunk<
  string,
  { runPostData: RunPostData },
  ThunkApiConfig
>(`${PIPELINE_SLICE_NAME}/run`, async ({ runPostData }, thunkAPI) => {
  const workspaceId = selectCurrentWorkspaceId(thunkAPI.getState())
  if (workspaceId) {
    try {
      const responseData = await runApi(workspaceId, runPostData)
      return responseData
    } catch (e) {
      return thunkAPI.rejectWithValue(e)
    }
  } else {
    return thunkAPI.rejectWithValue('workspace id does not exist.')
  }
})

export const runByCurrentUid = createAsyncThunk<
  string,
  { runPostData: Omit<RunPostData, 'name'> },
  ThunkApiConfig
>(
  `${PIPELINE_SLICE_NAME}/runByCurrentUid`,
  async ({ runPostData }, thunkAPI) => {
    const workspaceId = selectCurrentWorkspaceId(thunkAPI.getState())
    const currentUid = selectPipelineLatestUid(thunkAPI.getState())
    if (workspaceId && currentUid != null) {
      try {
        const responseData = await runByUidApi(
          workspaceId,
          currentUid,
          runPostData,
        )
        return responseData
      } catch (e) {
        return thunkAPI.rejectWithValue(e)
      }
    } else {
      return thunkAPI.rejectWithValue(
        'workspaceId or currentUid dose not exist.',
      )
    }
  },
)

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
  const workspaceId = selectCurrentWorkspaceId(thunkAPI.getState())
  if (workspaceId) {
    try {
      const responseData = await runResult({
        workspaceId,
        uid,
        pendingNodeIdList,
      })
      return responseData
    } catch (e) {
      return thunkAPI.rejectWithValue(e)
    }
  } else {
    return thunkAPI.rejectWithValue('workspace id does not exist')
  }
})
