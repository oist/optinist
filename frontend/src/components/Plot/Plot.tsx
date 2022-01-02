import React from 'react'
import { DisplayDataTabContext } from 'App'
import { DATA_TYPE_SET } from 'store/slice/DisplayData/DisplayDataType'
import { TablePlot } from './TablePlot'
import { TimeSeries } from './TimeSeries'
import { HeatMap } from './HeatMap'
import { ImagePlot } from './ImagePlot'

const Plot = React.memo(() => {
  const { dataType } = React.useContext(DisplayDataTabContext)
  switch (dataType) {
    case DATA_TYPE_SET.TABLE:
      return <TablePlot />
    case DATA_TYPE_SET.TIME_SERIES:
      return <TimeSeries />
    case DATA_TYPE_SET.HEAT_MAP:
      return <HeatMap />
    case DATA_TYPE_SET.IMAGE:
      return <ImagePlot />
  }
})

export default Plot
