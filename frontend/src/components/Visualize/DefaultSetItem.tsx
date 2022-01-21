import React from 'react'
import { useSelector } from 'react-redux'
import { DefaultPlot } from './Plot/DefaultPlot'
import { FilePathSelect } from './FilePathSelect'
import {
  selectDefaultSetFilePath,
  selectDefaultSetNodeId,
} from 'store/slice/VisualizeItem/VisualizeItemSelectors'

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
