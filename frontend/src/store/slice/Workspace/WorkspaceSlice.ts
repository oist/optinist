import { PayloadAction, createSlice } from '@reduxjs/toolkit'
import { WORKSPACE_SLICE_NAME, Workspace } from './WorkspaceType'
import { importExperimentByUid } from '../Experiments/ExperimentsActions'

const initialState: Workspace = {
  workspaces: [{ workspace_id: 'default' }],
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
    setCurrentWorkspace: (state, action: PayloadAction<string>) => {
      state.currentWorkspace.workspaceId = action.payload
    },
    clearCurrentWorkspace: (state) => {
      state.currentWorkspace = {
        selectedTab: 0,
      }
    },
  },
  extraReducers(builder) {
    builder.addCase(importExperimentByUid.fulfilled, (state, action) => {
      state.currentWorkspace.workspaceId = action.meta.arg.workspaceId
    })
    // TODO: add case for set loading on get workspaces pending
  },
})

export const { setCurrentWorkspace, clearCurrentWorkspace, setActiveTab } =
  workspaceSlice.actions
export default workspaceSlice.reducer
