import { RootState } from 'store/store'

export const selectWorkspace = (state: RootState) => state.workspace

export const selectCurrentWorkspaceId = (state: RootState) =>
  state.workspace.currentWorkspaceId

export const selectWorkspaceList = (state: RootState) =>
  state.workspace.workspaces

export const selectIsLoadingWorkspaceList = (state: RootState) =>
  state.workspace.loading
