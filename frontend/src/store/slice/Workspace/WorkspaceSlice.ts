import { PayloadAction, createSlice } from '@reduxjs/toolkit'
import { WORKSPACE_SLICE_NAME, Workspace } from './WorkspaceType'
import { reproduceWorkflow } from '../Workflow/WorkflowActions'

const initialState: Workspace = {
  workspaces: [{ workspace_id: 1 }],
  currentWorkspace: {
    selectedTab: 0,
  },
  loading: false,
}

export const workspaceSlice = createSlice({
  name: WORKSPACE_SLICE_NAME,
  initialState,
  reducers: {
    setActiveTab: (state, action: PayloadAction<number>) => {
      state.currentWorkspace.selectedTab = action.payload
    },
    setCurrentWorkspace: (state, action: PayloadAction<number>) => {
      state.currentWorkspace.workspaceId = action.payload
    },
    clearCurrentWorkspace: (state) => {
      state.currentWorkspace = {
        selectedTab: 0,
      }
    },
  },
  extraReducers(builder) {
    builder.addCase(reproduceWorkflow.fulfilled, (state, action) => {
      state.currentWorkspace.workspaceId = action.meta.arg.workspaceId
    })
    // TODO: add case for set loading on get workspaces pending
  },
})

export const { setCurrentWorkspace, clearCurrentWorkspace, setActiveTab } =
  workspaceSlice.actions
export default workspaceSlice.reducer
