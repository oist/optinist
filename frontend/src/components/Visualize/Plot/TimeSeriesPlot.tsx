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
  selectTimeSeriesStd,
  selectTimeSeriesXrange,
} from 'store/slice/DisplayData/DisplayDataSelectors'
import { getTimeSeriesDataById } from 'store/slice/DisplayData/DisplayDataActions'
import { TimeSeriesData } from 'store/slice/DisplayData/DisplayDataType'
import {
  selectTimeSeriesItemCheckedList,
  selectTimeSeriesItemDisplayNumbers,
  selectTimeSeriesItemOffset,
  selectTimeSeriesItemRefImageItemId,
  selectTimeSeriesItemShowGrid,
  selectTimeSeriesItemShowLine,
  selectTimeSeriesItemShowTickLabels,
  selectTimeSeriesItemSpan,
  selectTimeSeriesItemXrange,
  selectTimeSeriesItemZeroLine,
  selectVisualizeItemHeight,
  selectVisualizeItemWidth,
} from 'store/slice/VisualizeItem/VisualizeItemSelectors'
import { setTimeSeriesItemCheckedList } from 'store/slice/VisualizeItem/VisualizeItemSlice'
import createColormap from 'colormap'
import { setTimeSeriesItemDisplayNumbers } from 'store/slice/VisualizeItem/VisualizeItemSlice'
import { DisplayIndexMap } from 'store/slice/VisualizeItem/VisualizeItemType'

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
      dispatch(getTimeSeriesDataById({ path, index: 0 }))
    }
  }, [dispatch, isInitialized, path, displayNumbers])

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

  const dataXrange = useSelector(selectTimeSeriesXrange(path))
  const dataStd = useSelector(selectTimeSeriesStd(path))

  const offset = useSelector(selectTimeSeriesItemOffset(itemId))
  const span = useSelector(selectTimeSeriesItemSpan(itemId))
  const showgrid = useSelector(selectTimeSeriesItemShowGrid(itemId))
  const showline = useSelector(selectTimeSeriesItemShowLine(itemId))
  const showticklabels = useSelector(selectTimeSeriesItemShowTickLabels(itemId))
  const zeroline = useSelector(selectTimeSeriesItemZeroLine(itemId))
  const xrange = useSelector(selectTimeSeriesItemXrange(itemId))
  const refImageItemId = useSelector(selectTimeSeriesItemRefImageItemId(itemId))
  const displayNumbers = useSelector(
    selectTimeSeriesItemDisplayNumbers(itemId, refImageItemId),
  )
  const checkedList = useSelector(selectTimeSeriesItemCheckedList(itemId))
  const width = useSelector(selectVisualizeItemWidth(itemId))
  const height = useSelector(selectVisualizeItemHeight(itemId))

  const dataKeys = Object.keys(timeSeriesData)

  const colorScale = createColormap({
    colormap: 'jet',
    nshades: 100, //maxIndex >= 6 ? maxIndex : 6,
    format: 'hex',
    alpha: 1,
  })

  React.useEffect(() => {
    if (Object.keys(checkedList).length === 0 && dataKeys.length !== 0) {
      const checkedList: DisplayIndexMap = Object.fromEntries(
        dataKeys.map((v, i) => {
          if (i === 0) {
            return [v, true]
          }
          return [v, false]
        }),
      )
      dispatch(
        setTimeSeriesItemCheckedList({
          itemId,
          checkedList,
        }),
      )
    }
  }, [timeSeriesData, dispatch, itemId, checkedList, dataKeys])

  const data = React.useMemo(() => {
    if (timeSeriesData === null) {
      return {}
    }
    return Object.fromEntries(
      Object.keys(timeSeriesData).map((key, i) => {
        const num_key = parseInt(key)
        let y = Object.values(timeSeriesData[key])
        const new_i = Math.floor((i % 10) * 10 + i / 10) % 100
        const name = `(${String(parseInt(key) + 1)})`
        var error_y = {
          type: 'data',
          array: Object.keys(dataStd).includes(key)
            ? Object.values(dataStd[key])
            : null,
          visible: true,
        }
        if (offset) {
          error_y = {
            type: 'data',
            array: null,
            visible: true,
          }
        }
        var visible: string | boolean = 'legendonly'
        if (displayNumbers.includes(num_key)) {
          visible = true
        }
        if (displayNumbers.includes(num_key) && offset) {
          const activeIdx: number = displayNumbers.findIndex(
            (v) => v === num_key,
          )
          const mean: number = y.reduce((a, b) => a + b) / y.length
          const std: number =
            span *
            Math.sqrt(y.reduce((a, b) => a + Math.pow(b - mean, 2)) / y.length)
          y = y.map((value) => (value - mean) / (std + 1e-10) + activeIdx)
        }

        return [
          key,
          {
            name: name,
            x: dataXrange,
            y: y,
            visible: visible,
            line: { color: colorScale[new_i] },
            error_y: error_y,
          },
        ]
      }),
    )
  }, [
    timeSeriesData,
    displayNumbers,
    offset,
    span,
    colorScale,
    dataStd,
    dataXrange,
  ])

  const annotations = React.useMemo(() => {
    if (Object.keys(data).length !== 0) {
      return displayNumbers.map((value) => {
        return {
          x: Number(dataXrange[dataXrange.length - 1]) + dataXrange.length / 10,
          y: data[value.toString()].y[dataXrange.length - 1],
          xref: 'x',
          yref: 'y',
          text: `cell: ${value + 1}`,
          arrowhead: 1,
          ax: 0,
          ay: -10,
        }
      })
    } else {
      return []
    }
  }, [data, displayNumbers, dataXrange])

  const layout = React.useMemo(
    () => ({
      margin: {
        t: 60, // top
        l: 50, // left
        b: 30, // bottom
      },
      dragmode: 'pan',
      autosize: true,
      width: width,
      height: height - 50,
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
      annotations: annotations,
    }),
    [
      xrange,
      showgrid,
      showline,
      showticklabels,
      zeroline,
      annotations,
      width,
      height,
    ],
  )

  const config = {
    displayModeBar: true,
    responsive: true,
  }

  const onLegendClick = (event: LegendClickEvent) => {
    const clickNumber = parseInt(Object.keys(timeSeriesData)[event.curveNumber])

    const newDisplayNumbers = displayNumbers.includes(clickNumber)
      ? displayNumbers.filter((value) => value !== clickNumber)
      : [...displayNumbers, clickNumber]

    const newCheckedList = Object.fromEntries(
      Object.entries(checkedList).map(([key, value]) => {
        if (parseInt(key) === clickNumber) {
          if (displayNumbers.includes(clickNumber)) {
            return [key, false]
          } else {
            return [key, true]
          }
        }
        return [key, value]
      }),
    )

    dispatch(
      setTimeSeriesItemDisplayNumbers({
        itemId,
        displayNumbers: newDisplayNumbers,
      }),
    )

    dispatch(
      setTimeSeriesItemCheckedList({
        itemId,
        checkedList: newCheckedList,
      }),
    )

    // set DisplayNumbers
    if (!displayNumbers.includes(clickNumber)) {
      dispatch(getTimeSeriesDataById({ path, index: clickNumber }))
    }

    return false
  }

  return (
    <PlotlyChart
      data={Object.values(data)}
      layout={layout}
      config={config}
      onLegendClick={onLegendClick}
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
