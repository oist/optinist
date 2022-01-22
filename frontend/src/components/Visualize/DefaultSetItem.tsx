import React from 'react'
import { useSelector } from 'react-redux'
// import { DefaultPlot } from './Plot/DefaultPlot'
import { FilePathSelect } from './FilePathSelect'
import {
  selectDefaultSetFilePath,
  selectDefaultSetNodeId,
} from 'store/slice/VisualizeItem/VisualizeItemSelectors'
import { ImagePlot } from './Plot/ImagePlot'
import { DATA_TYPE } from 'store/slice/DisplayData/DisplayDataType'
import { DisplayDataContext } from './DataContext'

export const DefaultSetItem = React.memo<{
  itemId: number
}>(({ itemId }) => {
  return (
    <>
      <FilePathSelect
        itemId={itemId}
        dataType={'image'}
        selectedNodeId={useSelector(selectDefaultSetNodeId(itemId, 'image'))}
        selectedFilePath={useSelector(
          selectDefaultSetFilePath(itemId, 'image'),
        )}
      />
      <FilePathSelect
        itemId={itemId}
        dataType={'timeSeries'}
        selectedNodeId={useSelector(
          selectDefaultSetNodeId(itemId, 'timeSeries'),
        )}
        selectedFilePath={useSelector(
          selectDefaultSetFilePath(itemId, 'timeSeries'),
        )}
      />
      <FilePathSelect
        itemId={itemId}
        dataType={'heatMap'}
        selectedNodeId={useSelector(selectDefaultSetNodeId(itemId, 'heatMap'))}
        selectedFilePath={useSelector(
          selectDefaultSetFilePath(itemId, 'heatMap'),
        )}
      />
      <DefaultPlot itemId={itemId} />
    </>
  )
})

const DefaultPlot = React.memo<{
  itemId: number
}>(({ itemId }) => {
  const dataType = 'image'
  const filePath = useSelector(selectDefaultSetFilePath(itemId, dataType))
  const nodeId = useSelector(selectDefaultSetNodeId(itemId, dataType))
  if (filePath != null && dataType != null) {
    return (
      <DisplayDataContext.Provider
        value={{ nodeId, filePath, dataType, itemId }}
      >
        <ImagePlot />
        <div>default plot</div>
      </DisplayDataContext.Provider>
    )
  } else {
    return <div>Please select item correctly.</div>
  }
})
