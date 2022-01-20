import React, { useState } from 'react'
import { useSelector, useDispatch } from 'react-redux'
import PlotlyChart from 'react-plotlyjs-ts'
import { LinearProgress, Typography } from '@material-ui/core'

import { DisplayDataContext } from '../DisplayDataItem'
import {
  selectTimeSeriesData,
  selectTimeSeriesDataError,
  selectTimeSeriesDataIsFulfilled,
  selectTimeSeriesDataIsInitialized,
  selectTimeSeriesDataIsPending,
  selectTimeSeriesPlotlyData,
} from 'store/slice/DisplayData/DisplayDataSelectors'
import { getTimeSeriesData } from 'store/slice/DisplayData/DisplayDataActions'
import { TimeSeriesData } from 'store/slice/DisplayData/DisplayDataType'
import { selectNodeLabelById } from 'store/slice/FlowElement/FlowElementSelectors'
import { RootState } from 'store/store'
import { LegendClickEvent } from 'plotly.js'
import {
  selectTimeSeriesItemOffset,
  selectTimeSeriesItemShowGrid,
  selectTimeSeriesItemShowLine,
  selectTimeSeriesItemShowTickLabels,
  selectTimeSeriesItemSpan,
  selectTimeSeriesItemXrange,
  selectTimeSeriesItemZeroLine,
} from 'store/slice/VisualizeItem/VisualizeItemSelectors'

export const TimeSeries = React.memo(() => {
  const { filePath: path } = React.useContext(DisplayDataContext)
  const dispatch = useDispatch()
  const isPending = useSelector(selectTimeSeriesDataIsPending(path))
  const isInitialized = useSelector(selectTimeSeriesDataIsInitialized(path))
  const error = useSelector(selectTimeSeriesDataError(path))
  const isFulfilled = useSelector(selectTimeSeriesDataIsFulfilled(path))
  React.useEffect(() => {
    if (!isInitialized) {
      dispatch(getTimeSeriesData({ path }))
    }
  }, [dispatch, isInitialized, path])
  if (isPending) {
    return <LinearProgress />
  } else if (error != null) {
    return <Typography color="error">{error}</Typography>
  } else if (isFulfilled) {
    return <TimeSeriesImple />
  } else {
    return null
  }
})

const TimeSeriesImple = React.memo(() => {
  const { filePath: path, itemId } = React.useContext(DisplayDataContext)

  const timeSeriesData = useSelector(
    selectTimeSeriesData(path),
    timeSeriesDataEqualityFn,
  )

  const [displayNumbers, setDisplayNumbers] = useState([0])
  const offset = useSelector(selectTimeSeriesItemOffset(itemId))
  const span = useSelector(selectTimeSeriesItemSpan(itemId))
  const showgrid = useSelector(selectTimeSeriesItemShowGrid(itemId))
  const showline = useSelector(selectTimeSeriesItemShowLine(itemId))
  const showticklabels = useSelector(selectTimeSeriesItemShowTickLabels(itemId))
  const zeroline = useSelector(selectTimeSeriesItemZeroLine(itemId))
  const xrange = useSelector(selectTimeSeriesItemXrange(itemId))

  const data = React.useMemo(() => {
    if (timeSeriesData == null) {
      return []
    }
    return Object.keys(timeSeriesData['0'])
      .map((_, i) => {
        return {
          name: `(${i})`,
          x: Object.keys(timeSeriesData),
          y: Object.values(timeSeriesData).map((value) => value[i]),
        }
      })
      .map((v, idx) => {
        if (displayNumbers.includes(idx)) {
          if (offset) {
            const activeIdx: number = displayNumbers.findIndex(
              (value) => value == idx,
            )
            const mean: number = v.y.reduce((a, b) => a + b) / v.y.length
            const std: number =
              span *
              Math.sqrt(
                v.y.reduce((a, b) => a + Math.pow(b - mean, 2)) / v.y.length,
              )
            const newArray = v.y.map(
              (value) => (value - mean) / std + activeIdx,
            )
            return { name: v.name, x: v.x, y: newArray, visible: true }
          } else {
            return { name: v.name, x: v.x, y: v.y, visible: true }
          }
        } else {
          return { name: v.name, x: v.x, y: v.y, visible: 'legendonly' }
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
    [xrange, showgrid, showline, showticklabels, zeroline],
  )

  const config = {
    displayModeBar: true,
    scrollZoom: true,
  }

  const onClick = (event: any) => {
    const clickNumber = event.curveNumber
    if (displayNumbers.includes(clickNumber)) {
      setDisplayNumbers(displayNumbers.filter((value) => value != clickNumber))
    } else {
      setDisplayNumbers((prev) => [...prev, clickNumber])
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
