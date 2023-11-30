import { selectModeStandalone } from "store/slice/Standalone/StandaloneSeclector"
import { RootState } from "store/store"

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

export const selectRoiFilePathCancel = (state: RootState) =>
  state.workspace.currentWorkspace.roiFilePath

export const selectCurrentWorkspaceId = (state: RootState) =>
  state.workspace.currentWorkspace.workspaceId

export const selectStatusRoiCancel = (state: RootState) =>
  state.workspace.currentWorkspace.statusRoi

export const selectCurrentWorkspaceName = (state: RootState) =>
  state.workspace.currentWorkspace.workspaceName

export const selectCurrentWorkspaceOwnerId = (state: RootState) =>
  state.workspace.currentWorkspace.ownerId

export const selectIsLoadingWorkspaceList = (state: RootState) =>
  state.workspace.loading

export const selectIsWorkspaceOwner = (state: RootState) =>
  selectModeStandalone(state)
    ? true
    : state.workspace.currentWorkspace.ownerId === state.user.currentUser?.id
