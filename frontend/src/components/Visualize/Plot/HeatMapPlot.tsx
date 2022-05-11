import React from 'react'
import { useSelector, useDispatch } from 'react-redux'
import PlotlyChart from 'react-plotlyjs-ts'
import { LinearProgress, Typography } from '@mui/material'

import { DisplayDataContext } from '../DataContext'
import { twoDimarrayEqualityFn } from 'utils/EqualityUtils'
import {
  selectHeatMapData,
  selectHeatMapDataError,
  selectHeatMapDataIsFulfilled,
  selectHeatMapDataIsInitialized,
  selectHeatMapDataIsPending,
} from 'store/slice/DisplayData/DisplayDataSelectors'
import { getHeatMapData } from 'store/slice/DisplayData/DisplayDataActions'
import {
  selectHeatMapItemColors,
  selectHeatMapItemShowScale,
  selectVisualizeItemHeight,
  selectVisualizeItemWidth,
  selectVisualizeSaveFilename,
  selectVisualizeSaveFormat,
} from 'store/slice/VisualizeItem/VisualizeItemSelectors'

export const HeatMapPlot = React.memo(() => {
  const { filePath: path } = React.useContext(DisplayDataContext)
  const dispatch = useDispatch()
  const isPending = useSelector(selectHeatMapDataIsPending(path))
  const isInitialized = useSelector(selectHeatMapDataIsInitialized(path))
  const error = useSelector(selectHeatMapDataError(path))
  const isFulfilled = useSelector(selectHeatMapDataIsFulfilled(path))
  React.useEffect(() => {
    if (!isInitialized) {
      dispatch(getHeatMapData({ path }))
    }
  }, [dispatch, isInitialized, path])
  if (isPending) {
    return <LinearProgress />
  } else if (error != null) {
    return <Typography color="error">{error}</Typography>
  } else if (isFulfilled) {
    return <HeatMapImple />
  } else {
    return null
  }
})

const HeatMapImple = React.memo(() => {
  const { filePath: path, itemId } = React.useContext(DisplayDataContext)
  const heatMapData = useSelector(selectHeatMapData(path), heatMapDataEqualtyFn)
  const showscale = useSelector(selectHeatMapItemShowScale(itemId))
  const colorscale = useSelector(selectHeatMapItemColors(itemId))
  const width = useSelector(selectVisualizeItemWidth(itemId))
  const height = useSelector(selectVisualizeItemHeight(itemId))

  const data = React.useMemo(
    () =>
      heatMapData != null
        ? [
            {
              z: heatMapData,
              x: heatMapData.map((_, i) => i),
              y: heatMapData[0].map((_, i) => i),
              type: 'heatmap',
              name: 'heatmap',
              colorscale: colorscale.map((value) => {
                let offset: number = parseFloat(value.offset)
                const offsets: number[] = colorscale.map((v) => {
                  return parseFloat(v.offset)
                })
                // plotlyは端[0.0, 1.0]がないとダメなので、その設定
                if (offset === Math.max(...offsets)) {
                  offset = 1.0
                }
                if (offset === Math.min(...offsets)) {
                  offset = 0.0
                }
                return [offset, value.rgb]
              }),
              hoverongaps: false,
              showlegend: true,
              showscale: showscale,
            },
          ]
        : [],
    [heatMapData, showscale, colorscale],
  )

  const layout = React.useMemo(
    () => ({
      width: width,
      height: height - 50,
      dragmode: 'pan',
      margin: {
        t: 60, // top
        l: 50, // left
        b: 30, // bottom
      },
      autosize: true,
    }),
    [width, height],
  )

  const saveFileName = useSelector(selectVisualizeSaveFilename(itemId))
  const saveFormat = useSelector(selectVisualizeSaveFormat(itemId))

  const config = {
    displayModeBar: true,
    responsive: true,
    toImageButtonOptions: {
      format: saveFormat,
      filename: saveFileName,
    },
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
