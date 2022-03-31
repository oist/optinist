import React, { useState, useEffect } from 'react'

import { useSelector, useDispatch } from 'react-redux'
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
import { FilePathSelect } from './FilePathSelect'
import { setDisplayDataPath } from 'store/slice/VisualizeItem/VisualizeItemSlice'
import { deleteDisplayItem } from 'store/slice/DisplayData/DisplayDataSlice'

export const DisplayDataItem = React.memo<{
  itemId: number
}>(({ itemId }) => {
  const filePath = useSelector(selectVisualizeDataFilePath(itemId))
  const nodeId = useSelector(selectVisualizeDataNodeId(itemId))
  const dataType = useSelector(selectVisualizeDataType(itemId))

  const [prevItem, setPrevItem] = useState<{
    dataType: DATA_TYPE
    filePath: string | null
  }>({
    dataType: 'image',
    filePath: null,
  })

  useEffect(() => {
    setPrevItem({ dataType, filePath })
  }, [filePath, dataType])

  const dispatch = useDispatch()
  const onSelect = (nodeId: string, filePath: string, dataType: DATA_TYPE) => {
    dispatch(setDisplayDataPath({ itemId, nodeId, filePath, dataType }))
    dispatch(deleteDisplayItem(prevItem))
  }

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
