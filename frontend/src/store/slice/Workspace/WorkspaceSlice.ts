import { PayloadAction, createSlice } from '@reduxjs/toolkit'
import { WORKSPACE_SLICE_NAME, Workspace } from './WorkspaceType'
import { importExperimentByUid } from '../Experiments/ExperimentsActions'

const initialState: Workspace = {
  workspaces: [{ workspace_id: 'default' }],
  loading: false,
}

export const workspaceSlice = createSlice({
  name: WORKSPACE_SLICE_NAME,
  initialState,
  reducers: {
    setCurrentWorkspace: (state, action: PayloadAction<string>) => {
      state.currentWorkspaceId = action.payload
    },
    clearCurrentWorkspace: (state) => {
      state.currentWorkspaceId = undefined
    },
  },
  extraReducers(builder) {
    builder.addCase(importExperimentByUid.fulfilled, (state, action) => {
      state.currentWorkspaceId = action.meta.arg.workspaceId
    })
    // TODO: add case for set loading on get workspaces pending
  },
})

export const { setCurrentWorkspace, clearCurrentWorkspace } =
  workspaceSlice.actions
export default workspaceSlice.reducer
