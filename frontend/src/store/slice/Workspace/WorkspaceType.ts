export const WORKSPACE_SLICE_NAME = 'workspace'

export type Workspace = {
  workspaces: WorkspaceType[]
  currentWorkspace: {
    workspaceId?: string
    selectedTab: number
  }
  loading: boolean
}

export type WorkspaceType = {
  workspace_id: string
  // TODO: add fields required for workspace
}
