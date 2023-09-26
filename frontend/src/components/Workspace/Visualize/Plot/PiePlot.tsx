import React from 'react'
import { useSelector, useDispatch } from 'react-redux'
import PlotlyChart from 'react-plotlyjs-ts'
import { DisplayDataContext } from '../DataContext'
import {
  selectPieColumns,
  selectPieData,
  selectPieDataError,
  selectPieDataIsFulfilled,
  selectPieDataIsInitialized,
  selectPieDataIsPending,
} from 'store/slice/DisplayData/DisplayDataSelectors'
import { getPieData } from 'store/slice/DisplayData/DisplayDataActions'
import { LinearProgress, Typography } from '@mui/material'
import {
  selectVisualizeItemHeight,
  selectVisualizeItemWidth,
} from 'store/slice/VisualizeItem/VisualizeItemSelectors'

export const PiePlot = React.memo(() => {
  const { filePath: path } = React.useContext(DisplayDataContext)
  const dispatch = useDispatch()
  const isPending = useSelector(selectPieDataIsPending(path))
  const isInitialized = useSelector(selectPieDataIsInitialized(path))
  const error = useSelector(selectPieDataError(path))
  const isFulfilled = useSelector(selectPieDataIsFulfilled(path))

  React.useEffect(() => {
    if (!isInitialized) {
      dispatch(getPieData({ path }))
    }
  }, [dispatch, isInitialized, path])

  if (isPending) {
    return <LinearProgress />
  } else if (error != null) {
    return <Typography color="error">{error}</Typography>
  } else if (isFulfilled) {
    return <PiePlotImple />
  } else {
    return null
  }
})

const PiePlotImple = React.memo(() => {
  const { filePath: path, itemId } = React.useContext(DisplayDataContext)
  const pieData = useSelector(selectPieData(path))
  const columns = useSelector(selectPieColumns(path))
  const width = useSelector(selectVisualizeItemWidth(itemId))
  const height = useSelector(selectVisualizeItemHeight(itemId))

  const data = React.useMemo(
    () =>
      pieData != null
        ? [
            {
              values: pieData[0],
              labels: columns,
              type: 'pie',
              sort: false,
              direction: 'clockwise',
            },
          ]
        : [],
    [pieData, columns],
  )

  const layout = React.useMemo(
    () => ({
      width: width,
      height: height - 120,
      dragmode: 'pan',
      margin: {
        t: 60, // top
        l: 50, // left
        b: 30, // bottom
      },
      autosize: true,
      xaxis: {
        tickvals: columns,
      },
    }),
    [width, height, columns],
  )

  return <PlotlyChart data={data} layout={layout} />
})
