import { createContext } from "react"

import { DATA_TYPE } from "store/slice/DisplayData/DisplayDataType"

export const DisplayDataContext = createContext<{
  nodeId: string | null
  filePath: string
  dataType: DATA_TYPE
  itemId: number
}>({ nodeId: "", filePath: "", dataType: "csv", itemId: NaN })
