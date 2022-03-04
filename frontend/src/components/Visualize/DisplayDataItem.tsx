import React from 'react'

import { useSelector } from 'react-redux'
import {
  DATA_TYPE,
  DATA_TYPE_SET,
} from 'store/slice/DisplayData/DisplayDataType'
import {
  selectVisualizeDataFilePath,
  selectVisualizeDataNodeId,
  selectVisualizeDataType,
} from 'store/slice/VisualizeItem/VisualizeItemSelectors'
import { DisplayDataContext } from './DataContext'
import { HeatMapPlot } from './Plot/HeatMapPlot'
import { ImagePlot } from './Plot/ImagePlot'
import { RoiPlot } from './Plot/RoiPlot'
import { ScatterPlot } from './Plot/ScatterPlot'
import { CsvPlot } from './Plot/CsvPlot'
import { TimeSeriesPlot } from './Plot/TimeSeriesPlot'
import { BarPlot } from './Plot/BarPlot'
import { HTMLPlot } from './Plot/HTMLPlot'

export const DisplayDataItem = React.memo<{
  itemId: number
}>(({ itemId }) => {
  const filePath = useSelector(selectVisualizeDataFilePath(itemId))
  const nodeId = useSelector(selectVisualizeDataNodeId(itemId))
  const dataType = useSelector(selectVisualizeDataType(itemId))
  if (filePath != null && dataType != null) {
    return (
      <DisplayDataContext.Provider
        value={{ nodeId, filePath, dataType, itemId }}
      >
        <DisplayPlot dataType={dataType} />
      </DisplayDataContext.Provider>
    )
  } else {
    return <div>Please select item correctly.</div>
  }
})

const DisplayPlot = React.memo<{
  dataType: DATA_TYPE
}>(({ dataType }) => {
  switch (dataType) {
    case DATA_TYPE_SET.CSV:
      return <CsvPlot />
    case DATA_TYPE_SET.TIME_SERIES:
      return <TimeSeriesPlot />
    case DATA_TYPE_SET.HEAT_MAP:
      return <HeatMapPlot />
    case DATA_TYPE_SET.IMAGE:
      return <ImagePlot />
    case DATA_TYPE_SET.ROI:
      return <RoiPlot />
    case DATA_TYPE_SET.SCATTER:
      return <ScatterPlot />
    case DATA_TYPE_SET.BAR:
      return <BarPlot />
    case DATA_TYPE_SET.HTML:
      return <HTMLPlot />
    default:
      return null
  }
})
