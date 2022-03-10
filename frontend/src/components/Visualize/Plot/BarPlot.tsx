import React, { useEffect } from 'react'
import { useSelector, useDispatch } from 'react-redux'
import PlotlyChart from 'react-plotlyjs-ts'
import { LinearProgress, Typography } from '@mui/material'

import { DisplayDataContext } from '../DataContext'
import {
  selectBarData,
  selectBarDataError,
  selectBarDataIsFulfilled,
  selectBarDataIsInitialized,
  selectBarDataIsPending,
} from 'store/slice/DisplayData/DisplayDataSelectors'
import { getBarData } from 'store/slice/DisplayData/DisplayDataActions'
import { BarData } from 'store/slice/DisplayData/DisplayDataType'

export const BarPlot = React.memo(() => {
  const { filePath: path } = React.useContext(DisplayDataContext)
  const dispatch = useDispatch()
  const isPending = useSelector(selectBarDataIsPending(path))
  const isInitialized = useSelector(selectBarDataIsInitialized(path))
  const error = useSelector(selectBarDataError(path))
  const isFulfilled = useSelector(selectBarDataIsFulfilled(path))
  React.useEffect(() => {
    if (!isInitialized) {
      dispatch(getBarData({ path }))
    }
  }, [dispatch, isInitialized, path])
  if (isPending) {
    return <LinearProgress />
  } else if (error != null) {
    return <Typography color="error">{error}</Typography>
  } else if (isFulfilled) {
    return <BarPlotImple />
  } else {
    return null
  }
})

const BarPlotImple = React.memo(() => {
  const { filePath: path } = React.useContext(DisplayDataContext)

  const barData = useSelector(selectBarData(path), barDataEqualityFn)

  const data = React.useMemo(
    () => [
      {
        x: Object.keys(barData[0]),
        y: Object.values(barData[0]),
        type: 'bar',
      },
    ],
    [barData],
  )

  const layout = React.useMemo(
    () => ({
      title: path.split('/').reverse()[0],
      margin: {
        t: 60, // top
        l: 50, // left
        b: 30, // bottom
      },
      dragmode: 'pan',
      autosize: true,
    }),
    [path],
  )

  const config = {
    displayModeBar: true,
    // scrollZoom: true,
    responsive: true,
  }

  const ref = React.useRef<HTMLDivElement>(null)
  const plotlyHeight = ref.current?.getBoundingClientRect().height

  useEffect(() => {
    const height =
      ref.current?.getElementsByClassName('main-svg')[0].clientHeight
    const plotContainer = (
      ref.current?.getElementsByClassName(
        'plot-container',
      ) as HTMLCollectionOf<HTMLElement>
    )[0]

    if (height !== undefined && plotContainer !== undefined) {
      plotContainer.style.height = `${height}px`
    }
  }, [plotlyHeight])

  return (
    <div ref={ref}>
      <PlotlyChart data={data} layout={layout} config={config} />
    </div>
  )
})

function barDataEqualityFn(a: BarData | undefined, b: BarData | undefined) {
  if (a != null && b != null) {
    const aArray = Object.keys(a)
    const bArray = Object.keys(b)
    return (
      a === b ||
      (aArray.length === bArray.length &&
        aArray.every((aKey, i) => {
          const bKey = bArray[i]
          return bKey === aKey // && nestEqualityFun(bValue, aValue)
        }))
    )
  } else {
    return a === undefined && b === undefined
  }
}
