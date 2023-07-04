import { createSlice } from '@reduxjs/toolkit'
import { WORKSPACE_SLICE_NAME, Workspace } from './WorkspaceType'
import { importExperimentByUid } from '../Experiments/ExperimentsActions'

const initialState: Workspace = {
  workspaces: [],
}

export const workspaceSlice = createSlice({
  name: WORKSPACE_SLICE_NAME,
  initialState,
  reducers: {},
  extraReducers(builder) {
    builder.addCase(importExperimentByUid.fulfilled, (state, action) => {
      state.currentWorkspace = {
        workspace_id: action.meta.arg.workspaceId,
      }
    })
  },
})

export default workspaceSlice.reducer
