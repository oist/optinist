import { RunPostData } from 'api/run/Run'
import { WORKFLOW_SLICE_NAME } from './WorkflowType'
import { createAsyncThunk } from '@reduxjs/toolkit'
import {
  importWorkflowByUidApi,
  loadWorkflowConfigApi,
} from 'api/workflow/Workflow'

export const importWorkflowByUid = createAsyncThunk<
  RunPostData,
  { workspaceId: string; uid: string }
>(
  `${WORKFLOW_SLICE_NAME}/importWorkflowByUid`,
  async ({ workspaceId, uid }, thunkAPI) => {
    try {
      const response = await importWorkflowByUidApi(workspaceId, uid)
      return response
    } catch (e) {
      return thunkAPI.rejectWithValue(e)
    }
  },
)

export const loadWorkflowConfig = createAsyncThunk<
  RunPostData,
  { formData: FormData }
>(
  `${WORKFLOW_SLICE_NAME}/loadWorkflowConfig`,
  async ({ formData }, thunkAPI) => {
    try {
      const response = await loadWorkflowConfigApi(formData)
      return response
    } catch (e) {
      return thunkAPI.rejectWithValue(e)
    }
  },
)
