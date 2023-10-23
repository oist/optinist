import { createAsyncThunk } from "@reduxjs/toolkit"

import {
  fetchWorkflowApi,
  reproduceWorkflowApi,
  importWorkflowConfigApi,
  WorkflowConfigDTO,
  WorkflowWithResultDTO,
} from "api/workflow/Workflow"
import { WORKFLOW_SLICE_NAME } from "store/slice/Workflow/WorkflowType"

export const fetchWorkflow = createAsyncThunk<WorkflowWithResultDTO, number>(
  `${WORKFLOW_SLICE_NAME}/fetchExperiment`,
  async (workspaceId, thunkAPI) => {
    try {
      const response = await fetchWorkflowApi(workspaceId)
      return response
    } catch (e) {
      return thunkAPI.rejectWithValue(e)
    }
  },
)

export const reproduceWorkflow = createAsyncThunk<
  WorkflowWithResultDTO,
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
