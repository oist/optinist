export enum WORKSPACE_TYPE {
  DEFAULT = 0,
  NORMAL = 1,
}

export const WORKSPACE_TYPE_LABEL: Record<WORKSPACE_TYPE, string> = {
  [WORKSPACE_TYPE.DEFAULT]: "Normal", // Same as NORMAL
  [WORKSPACE_TYPE.NORMAL]: "Normal",
}
