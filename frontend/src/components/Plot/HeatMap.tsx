import React from 'react'

import { useSelector, useDispatch } from 'react-redux'

import PlotlyChart from 'react-plotlyjs-ts'

import { OutputPlotContext } from 'App'
import {
  heatMapDataSelector,
  heatMapDataIsLoadedSelector,
} from 'store/slice/PlotData/PlotDataSelector'
import { outputPathValueByIdSelector } from 'store/slice/Algorithm/AlgorithmSelector'
import { getHeatMapData } from 'store/slice/PlotData/PlotDataAction'
import { twoDimarrayEqualityFn } from 'utils/EqualityUtils'

export const HeatMap = React.memo(() => {
  const { nodeId, outputKey } = React.useContext(OutputPlotContext)
  const dispatch = useDispatch()
  const path = useSelector(outputPathValueByIdSelector(nodeId, outputKey))
  const isLoaded = useSelector(heatMapDataIsLoadedSelector(path ?? ''))
  React.useEffect(() => {
    if (!isLoaded && path != null) {
      dispatch(getHeatMapData({ path }))
    }
  }, [dispatch, isLoaded, path])
  if (isLoaded) {
    return <HeatMapImple />
  } else {
    return null
  }
})
const HeatMapImple = React.memo(() => {
  const { nodeId, outputKey } = React.useContext(OutputPlotContext)
  const path = useSelector(outputPathValueByIdSelector(nodeId, outputKey))
  const heatMapData = useSelector(
    heatMapDataSelector(path ?? ''),
    heatMapDataEqualtyFn,
  )

  const data = React.useMemo(
    () =>
      heatMapData != null
        ? [
            {
              z: heatMapData,
              x: heatMapData.map((_, i) => i),
              y: heatMapData[0].map((_, i) => i),
              type: 'heatmap',
              hoverongaps: false,
            },
          ]
        : [],
    [heatMapData],
  )

  const layout = {
    title: 'dummy',
    margin: {
      t: 60, // top
      l: 50, // left
      b: 30, // bottom
    },
    autosize: true,
    height: 350,
  }
  const config = {
    displayModeBar: true,
  }
  return <PlotlyChart data={data} layout={layout} config={config} />
})

function heatMapDataEqualtyFn(
  a: number[][] | undefined,
  b: number[][] | undefined,
) {
  if (a != null && b != null) {
    return twoDimarrayEqualityFn(a, b)
  } else {
    return a === undefined && b === undefined
  }
}
