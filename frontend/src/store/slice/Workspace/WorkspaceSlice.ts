import { PayloadAction, createSlice, isAnyOf } from "@reduxjs/toolkit"

import { reproduceWorkflow } from "store/slice/Workflow/WorkflowActions"
import {
  delWorkspace,
  getListUserShareWorkSpaces,
  getWorkspace,
  getWorkspaceList,
  postListUserShareWorkspaces,
  postWorkspace,
  putWorkspace,
} from "store/slice/Workspace/WorkspaceActions"
import {
  WORKSPACE_SLICE_NAME,
  Workspace,
} from "store/slice/Workspace/WorkspaceType"

const initialState: Workspace = {
  currentWorkspace: {
    selectedTab: 0,
  },
  workspace: {
    items: [],
    total: 0,
    limit: 50,
    offset: 0,
  },
  loading: false,
  listUserShare: undefined,
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
    builder
      .addCase(reproduceWorkflow.fulfilled, (state, action) => {
        state.currentWorkspace.workspaceId = action.meta.arg.workspaceId
        state.currentWorkspace.selectedTab = 0
      })
      .addCase(getWorkspace.fulfilled, (state, action) => {
        state.currentWorkspace.workspaceId = action.payload.id
        state.currentWorkspace.workspaceName = action.payload.name
        state.currentWorkspace.ownerId = action.payload.user.id
        state.loading = false
      })
      .addCase(getWorkspaceList.fulfilled, (state, action) => {
        state.workspace = action.payload
        state.loading = false
      })
      .addCase(getListUserShareWorkSpaces.fulfilled, (state, action) => {
        state.listUserShare = action.payload
        state.loading = false
      })
      .addMatcher(
        isAnyOf(
          getWorkspace.rejected,
          getWorkspaceList.rejected,
          postWorkspace.fulfilled,
          postWorkspace.rejected,
          putWorkspace.fulfilled,
          putWorkspace.rejected,
          delWorkspace.fulfilled,
          delWorkspace.rejected,
          getListUserShareWorkSpaces.rejected,
          postListUserShareWorkspaces.rejected,
        ),
        (state) => {
          state.loading = false
        },
      )
      .addMatcher(
        isAnyOf(
          getWorkspace.pending,
          getWorkspaceList.pending,
          postWorkspace.pending,
          putWorkspace.pending,
          delWorkspace.pending,
          getListUserShareWorkSpaces.pending,
          postListUserShareWorkspaces.pending,
        ),
        (state) => {
          state.loading = true
        },
      )
  },
})

export const { setCurrentWorkspace, clearCurrentWorkspace, setActiveTab } =
  workspaceSlice.actions
export default workspaceSlice.reducer
