import React from 'react'
import { DATA_TYPE } from 'store/slice/DisplayData/DisplayDataType'

export const DisplayDataContext = React.createContext<{
  nodeId: string | null
  filePath: string
  dataType: DATA_TYPE
  itemId: number
}>({ nodeId: '', filePath: '', dataType: 'table', itemId: NaN })

export const RunPipeLineContext = React.createContext<{
  // const [triggerRunPipeline, result] = useLazyRunPipelineQuery()
  runPipeLine: any
  result: any
}>({ runPipeLine: null, result: null })
