import React from 'react'
import { useSelector } from 'react-redux'
import { DATA_TYPE } from 'store/slice/DisplayData/DisplayDataType'
import {
  selectDefaultSetFilePath,
  selectDefaultSetNodeId,
} from 'store/slice/VisualizeItem/VisualizeItemSelectors'
import { ImagePlot } from './ImagePlot'

export const DefaultPlot = React.memo<{
  itemId: number
}>(({ itemId }) => {
  const dataType = 'image'
  const filePath = useSelector(selectDefaultSetFilePath(itemId, dataType))
  const nodeId = useSelector(selectDefaultSetNodeId(itemId, dataType))
  return (
    <DisplayDataContext.Provider value={{ nodeId, filePath, dataType, itemId }}>
      <ImagePlot />
      <div>default plot</div>
    </DisplayDataContext.Provider>
  )
})

export const DisplayDataContext = React.createContext<{
  nodeId: string | null
  filePath: string | null
  dataType: DATA_TYPE | null
  itemId: number
}>({ nodeId: '', filePath: '', dataType: 'table', itemId: NaN })
