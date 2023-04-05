import { createContext } from 'react'

export const DialogContext = createContext<{
  onOpen: (nodeId: string) => any
}>({
  onOpen: () => null,
})
