import { IS_STANDALONE } from 'const/Mode'
import { RootState } from 'store/store'

export const selectWorkspace = (state: RootState) => state.workspace
export const selectWorkspaceListUserShare = (state: RootState) =>
  state.workspace.listUserShare
export const selectWorkspaceData = (state: RootState) =>
  state.workspace.workspace

export const selectWorkspaceItem =
  (workspaceId: number) => (state: RootState) =>
    selectWorkspaceData(state).items.find((item) => item.id === workspaceId)

export const selectActiveTab = (state: RootState) =>
  state.workspace.currentWorkspace.selectedTab

export const selectCurrentWorkspaceId = (state: RootState) =>
  state.workspace.currentWorkspace.workspaceId

export const selectCurrentWorkspaceOwnerId = (state: RootState) =>
  state.workspace.currentWorkspace.ownerId

export const selectIsLoadingWorkspaceList = (state: RootState) =>
  state.workspace.loading

export const selectIsWorkspaceOwner = (state: RootState) =>
  IS_STANDALONE
    ? true
    : state.workspace.currentWorkspace.ownerId === state.user.currentUser?.id
