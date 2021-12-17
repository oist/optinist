import React from 'react'
import { useSelector } from 'react-redux'

import { outputPathTypeByIdSelector } from 'store/slice/Algorithm/AlgorithmSelector'
import { OutputPlotContext } from 'App'
import { TimeSeries } from './Plot/TimeSeries'
import { HeatMap } from './Plot/HeatMap'
import { OUTPUT_TYPE_SET } from 'store/slice/Algorithm/AlgorithmType'
import 'style/plot.css'

export const Plot = React.memo(function PlotOutput() {
  const { nodeId, outputKey } = React.useContext(OutputPlotContext)
  const outputType = useSelector(outputPathTypeByIdSelector(nodeId, outputKey))
  switch (outputType) {
    case OUTPUT_TYPE_SET.HEAT_MAP:
      return <HeatMap />
    case OUTPUT_TYPE_SET.TIME_SERIES:
      return <TimeSeries />
    case OUTPUT_TYPE_SET.IMAGE:
    default:
      return null
  }
})
