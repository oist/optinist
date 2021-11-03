import React from 'react'
import { useSelector, useDispatch } from 'react-redux'

import { getAlgoOutputData } from 'redux/slice/Algorithm/AlgorithmAction'
import {
  outputDataIsLoadedByIdSelector,
  outputPathTypeByIdSelector,
  outputPathValueByIdSelector,
} from 'redux/slice/Algorithm/AlgorithmSelector'
import { OutputPlotContext } from 'App'
import { toOutputDataId } from 'redux/slice/Algorithm/AlgorithmUtils'
import { TimeSeries } from './TimeSeries'
import { HeatMap } from './HeatMap'
import { OUTPUT_TYPE_SET } from 'redux/slice/Algorithm/AlgorithmType'

export const Plot = React.memo(function PlotOutput() {
  const { nodeId, outputKey } = React.useContext(OutputPlotContext)
  const dispatch = useDispatch()
  const path = useSelector(outputPathValueByIdSelector(nodeId, outputKey))
  const isLoaded = useSelector(
    outputDataIsLoadedByIdSelector(toOutputDataId(nodeId, outputKey)),
  )
  React.useEffect(() => {
    if (!isLoaded && path != null) {
      dispatch(getAlgoOutputData({ nodeId, outputKey, path }))
    }
  }, [isLoaded, nodeId, path])
  const outputType = useSelector(outputPathTypeByIdSelector(nodeId, outputKey))
  if (isLoaded) {
    switch (outputType) {
      case OUTPUT_TYPE_SET.HEAT_MAP:
        return <HeatMap />
      case OUTPUT_TYPE_SET.TIME_SERIES:
        return <TimeSeries />
      case OUTPUT_TYPE_SET.IMAGE:
      default:
        return null
    }
  } else {
    return null
  }
})
