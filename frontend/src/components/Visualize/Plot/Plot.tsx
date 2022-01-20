import React from 'react'
import { DATA_TYPE_SET } from 'store/slice/DisplayData/DisplayDataType'
import { TablePlot } from './TablePlot'
import { TimeSeries } from './TimeSeries'
import { HeatMap } from './HeatMap'
import { ImagePlot } from './ImagePlot'
import { DisplayDataContext } from '../DisplayDataItem'
import { RoiPlot } from './RoiPlot'

const Plot = React.memo(() => {
  const { dataType } = React.useContext(DisplayDataContext)
  switch (dataType) {
    case DATA_TYPE_SET.TABLE:
      return <TablePlot />
    case DATA_TYPE_SET.TIME_SERIES:
      return <TimeSeries />
    case DATA_TYPE_SET.HEAT_MAP:
      return <HeatMap />
    case DATA_TYPE_SET.IMAGE:
      return <ImagePlot />
    case DATA_TYPE_SET.ROI:
      return <RoiPlot />
  }
})

export default Plot
