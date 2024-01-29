import { createContext } from "react"

import { FILE_TREE_TYPE } from "api/files/Files"

export declare type FileSelectDialogValue = {
  filePath: string | string[]
  open: boolean
  fileTreeType?: FILE_TREE_TYPE
  multiSelect: boolean
  onSelectFile: (path: string | string[]) => void
}

export declare type FileInputUrl = {
  filePath?: string | string[]
  open: boolean
  nodeId: string
  requestId: string
  fileType?: FILE_TREE_TYPE
}

export declare type ClearWorkflowIdDialogValue = {
  open: boolean
  handleOk: () => void
  handleCancel: () => void
}

export declare type ErrorDialogValue = {
  anchorElRef: { current: Element | null }
  message: string
}

export const DialogContext = createContext<{
  onOpenOutputDialog: (nodeId: string) => void
  onOpenFileSelectDialog: (value: FileSelectDialogValue) => void
  onOpenInputUrlDialog: (value: FileInputUrl) => void
  onOpenClearWorkflowIdDialog: (value: ClearWorkflowIdDialogValue) => void
  onMessageError: (value: ErrorDialogValue) => void
}>({
  onOpenOutputDialog: () => null,
  onOpenFileSelectDialog: () => null,
  onOpenClearWorkflowIdDialog: () => null,
  onOpenInputUrlDialog: () => null,
  onMessageError: () => null,
})
