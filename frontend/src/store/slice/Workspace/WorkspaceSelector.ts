import { RootState } from 'store/store'

export const selectWorkspace = (state: RootState) => state.workspace
// NOTE: optinist is still a single user app, so just has 'default' workspace.
export const selectCurrentWorkspaceId = (state: RootState) =>
  state.workspace.currentWorkspace?.workspace_id ?? 'default'
export const selectWorkspaceList = (state: RootState) =>
  state.workspace.workspaces
