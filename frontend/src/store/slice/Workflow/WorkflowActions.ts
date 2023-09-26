import { WORKFLOW_SLICE_NAME } from './WorkflowType'
import { createAsyncThunk } from '@reduxjs/toolkit'
import {
  reproduceWorkflowApi,
  importWorkflowConfigApi,
  WorkflowConfigDTO,
} from 'api/workflow/Workflow'

export const reproduceWorkflow = createAsyncThunk<
  WorkflowConfigDTO,
  { workspaceId: number; uid: string }
>(
  `${WORKFLOW_SLICE_NAME}/reproduceWorkflow`,
  async ({ workspaceId, uid }, thunkAPI) => {
    try {
      const response = await reproduceWorkflowApi(workspaceId, uid)
      return response
    } catch (e) {
      return thunkAPI.rejectWithValue(e)
    }
  },
)

export const importWorkflowConfig = createAsyncThunk<
  WorkflowConfigDTO,
  { formData: FormData }
>(
  `${WORKFLOW_SLICE_NAME}/importWorkflowConfig`,
  async ({ formData }, thunkAPI) => {
    try {
      const response = await importWorkflowConfigApi(formData)
      return response
    } catch (e) {
      return thunkAPI.rejectWithValue(e)
    }
  },
)
