import { createContext } from 'react'

export const DialogContext = createContext<{
  onOpen: (nodeId: string) => any
  onOpenDialogFile: (value: {
    filePath: string | string[]
    open: boolean
    fileTreeType?: string
    multiSelect: boolean
    onSelectFile: (v: any) => any
  }) => any
  onMessageError: (v: { anchorElRef: { current: Element | null }; message: string }) => any
}>({
  onOpen: () => null,
  onOpenDialogFile: () => null,
  onMessageError: () => null,
})
