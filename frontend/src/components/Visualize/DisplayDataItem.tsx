import React from 'react'

import { useSelector, useDispatch } from 'react-redux'
import { DATA_TYPE } from 'store/slice/DisplayData/DisplayDataType'
import {
  selectVisualizeDataFilePath,
  selectVisualizeDataNodeId,
  selectVisualizeDataType,
} from 'store/slice/VisualizeItem/VisualizeItemSelectors'
import { Plot } from './Plot'

export const DisplayDataItem = React.memo<{ itemId: number }>(({ itemId }) => {
  const filePath = useSelector(selectVisualizeDataFilePath(itemId))
  const nodeId = useSelector(selectVisualizeDataNodeId(itemId))
  const dataType = useSelector(selectVisualizeDataType(itemId))
  if (nodeId != null && filePath != null && dataType != null) {
    return (
      <DisplayDataContext.Provider value={{ nodeId, filePath, dataType }}>
        <Plot />
      </DisplayDataContext.Provider>
    )
  } else {
    return <div>Please select item correctly.</div>
  }
})

export const DisplayDataContext = React.createContext<{
  nodeId: string
  filePath: string
  dataType: DATA_TYPE
}>({ nodeId: '', filePath: '', dataType: 'table' })
