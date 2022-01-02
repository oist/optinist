import React from 'react'
import { useSelector, useDispatch } from 'react-redux'
import PlotlyChart from 'react-plotlyjs-ts'
import { LinearProgress, Typography } from '@material-ui/core'

import { DisplayDataTabContext } from 'App'
import { twoDimarrayEqualityFn } from 'utils/EqualityUtils'
import {
  selectHeatMapData,
  selectHeatMapDataError,
  selectHeatMapDataIsFulfilled,
  selectHeatMapDataIsInitialized,
  selectHeatMapDataIsPending,
} from 'store/slice/DisplayData/DisplayDataSelectors'
import { getHeatMapData } from 'store/slice/DisplayData/DisplayDataActions'
import { selectNodeLabelById } from 'store/slice/FlowElement/FlowElementSelectors'

export const HeatMap = React.memo(() => {
  const { filePath: path } = React.useContext(DisplayDataTabContext)
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
  const { filePath: path, nodeId } = React.useContext(DisplayDataTabContext)
  const label = useSelector(selectNodeLabelById(nodeId))
  const heatMapData = useSelector(selectHeatMapData(path), heatMapDataEqualtyFn)
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
    title: label,
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
