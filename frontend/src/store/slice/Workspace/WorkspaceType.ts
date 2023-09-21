export const WORKSPACE_SLICE_NAME = 'workspace'

export type Workspace = {
  workspaces: WorkspaceType[]
  currentWorkspace: {
    workspaceId?: number
    selectedTab: number
  }
  loading: boolean
}

export type WorkspaceType = {
  workspace_id: number
  // TODO: add fields required for workspace
}
