import React from 'react'
import { useSelector, useDispatch } from 'react-redux'
import PlotlyChart from 'react-plotlyjs-ts'
import { LegendClickEvent } from 'plotly.js'
import { LinearProgress, Typography } from '@mui/material'

import { DisplayDataContext } from '../DataContext'
import {
  selectTimeSeriesData,
  selectTimeSeriesDataError,
  selectTimeSeriesDataIsFulfilled,
  selectTimeSeriesDataIsInitialized,
  selectTimeSeriesDataIsPending,
} from 'store/slice/DisplayData/DisplayDataSelectors'
import { getTimeSeriesData } from 'store/slice/DisplayData/DisplayDataActions'
import { TimeSeriesData } from 'store/slice/DisplayData/DisplayDataType'
import {
  selectTimeSeriesItemDisplayNumbers,
  selectTimeSeriesItemOffset,
  selectTimeSeriesItemShowGrid,
  selectTimeSeriesItemShowLine,
  selectTimeSeriesItemShowTickLabels,
  selectTimeSeriesItemSpan,
  selectTimeSeriesItemXrange,
  selectTimeSeriesItemZeroLine,
} from 'store/slice/VisualizeItem/VisualizeItemSelectors'
import createColormap from 'colormap'
import { setTimeSeriesItemDisplayNumbers } from 'store/slice/VisualizeItem/VisualizeItemSlice'

export const TimeSeriesPlot = React.memo(() => {
  const { itemId, filePath: path } = React.useContext(DisplayDataContext)
  const dispatch = useDispatch()
  const displayNumbers = useSelector(selectTimeSeriesItemDisplayNumbers(itemId))
  const isPending = useSelector(selectTimeSeriesDataIsPending(path))
  const isInitialized = useSelector(selectTimeSeriesDataIsInitialized(path))
  const error = useSelector(selectTimeSeriesDataError(path))
  const isFulfilled = useSelector(selectTimeSeriesDataIsFulfilled(path))
  React.useEffect(() => {
    if (!isInitialized) {
      displayNumbers.map((v) => dispatch(getTimeSeriesData({ path, index: v })))
    }
  }, [dispatch, isInitialized, path])

  if (!isInitialized) {
    return <LinearProgress />
  } else if (error != null) {
    return <Typography color="error">{error}</Typography>
  } else if (isPending || isFulfilled) {
    return <TimeSeriesPlotImple />
  } else {
    return null
  }
})

const TimeSeriesPlotImple = React.memo(() => {
  const { filePath: path, itemId } = React.useContext(DisplayDataContext)

  // 0番のデータとkeysだけをとってくる
  const dispatch = useDispatch()
  const timeSeriesData = useSelector(
    selectTimeSeriesData(path),
    timeSeriesDataEqualityFn,
  )

  const offset = useSelector(selectTimeSeriesItemOffset(itemId))
  const span = useSelector(selectTimeSeriesItemSpan(itemId))
  const showgrid = useSelector(selectTimeSeriesItemShowGrid(itemId))
  const showline = useSelector(selectTimeSeriesItemShowLine(itemId))
  const showticklabels = useSelector(selectTimeSeriesItemShowTickLabels(itemId))
  const zeroline = useSelector(selectTimeSeriesItemZeroLine(itemId))
  const xrange = useSelector(selectTimeSeriesItemXrange(itemId))
  const displayNumbers = useSelector(selectTimeSeriesItemDisplayNumbers(itemId))

  const dataSize = Object.keys(timeSeriesData).length

  const colorScale = createColormap({
    colormap: 'jet',
    nshades: dataSize < 6 ? 6 : dataSize,
    format: 'hex',
    alpha: 1,
  })

  const data = React.useMemo(() => {
    if (timeSeriesData === null) {
      return []
    }
    return Object.keys(timeSeriesData).map((key, i) => {
      let y = Object.values(timeSeriesData[key])
      const new_i = (i * 10 + i) % dataSize

      if (displayNumbers.includes(i)) {
        if (offset) {
          const activeIdx: number = displayNumbers.findIndex(
            (value) => value === i,
          )
          const mean: number = y.reduce((a, b) => a + b) / y.length
          const std: number =
            span *
            Math.sqrt(y.reduce((a, b) => a + Math.pow(b - mean, 2)) / y.length)
          const newArray = y.map((value) => (value - mean) / std + activeIdx)
          return {
            name: `(${key})`,
            x: Object.keys(timeSeriesData[key]),
            y: newArray,
            visible: true,
            line: { color: colorScale[new_i] },
          }
        } else {
          return {
            name: `(${key})`,
            x: Object.keys(timeSeriesData[key]),
            y: y,
            visible: true,
            line: { color: colorScale[new_i] },
          }
        }
      } else {
        return {
          name: `(${key})`,
          x: Object.keys(timeSeriesData[key]),
          y: y,
          visible: 'legendonly',
          line: { color: colorScale[new_i] },
        }
      }
    })
  }, [timeSeriesData, displayNumbers, offset, span])

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
      height: 300,
      xaxis: {
        range: [xrange.left, xrange.right],
        showgrid: showgrid,
        showline: showline,
        showticklabels: showticklabels,
        zeroline: zeroline,
      },
      yaxis: {
        showgrid: showgrid,
        showline: showline,
        showticklabels: showticklabels,
        zeroline: zeroline,
      },
    }),
    [path, xrange, showgrid, showline, showticklabels, zeroline],
  )

  const config = {
    displayModeBar: true,
    responsive: true,
  }

  const onClick = (event: LegendClickEvent) => {
    const clickNumber = event.curveNumber
    if (displayNumbers.includes(clickNumber)) {
      const newValue = displayNumbers.filter((value) => value !== clickNumber)
      dispatch(
        setTimeSeriesItemDisplayNumbers({
          itemId,
          displayNumbers: newValue,
        }),
      )
    } else {
      const newValue = [...displayNumbers, clickNumber]
      dispatch(
        setTimeSeriesItemDisplayNumbers({
          itemId,
          displayNumbers: newValue,
        }),
      )
      dispatch(getTimeSeriesData({ path, index: clickNumber }))
    }
    return false
  }

  return (
    <PlotlyChart
      data={data}
      layout={layout}
      config={config}
      onLegendClick={onClick}
    />
  )
})

function timeSeriesDataEqualityFn(
  a: TimeSeriesData | undefined,
  b: TimeSeriesData | undefined,
) {
  if (a != null && b != null) {
    const aArray = Object.entries(a)
    const bArray = Object.entries(b)
    return (
      a === b ||
      (aArray.length === bArray.length &&
        aArray.every(([aKey, aValue], i) => {
          const [bKey, bValue] = bArray[i]
          return bKey === aKey && nestEqualityFun(bValue, aValue)
        }))
    )
  } else {
    return a === undefined && b === undefined
  }
}

function nestEqualityFun(
  a: {
    [key: number]: number
  },
  b: {
    [key: number]: number
  },
) {
  const aArray = Object.entries(a)
  const bArray = Object.entries(b)
  return (
    a === b ||
    (aArray.length === bArray.length &&
      aArray.every(([aKey, aValue], i) => {
        const [bKey, bValue] = bArray[i]
        return bKey === aKey && bValue === aValue
      }))
  )
}
