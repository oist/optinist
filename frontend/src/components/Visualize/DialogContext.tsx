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
}>({
  onOpen: () => null,
  onOpenDialogFile: () => null,
})
