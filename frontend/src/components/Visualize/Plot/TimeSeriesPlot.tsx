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
import {
  getTimeSeriesDataById,
  getTimeSeriesInitData,
} from 'store/slice/DisplayData/DisplayDataActions'
import { TimeSeriesData } from 'api/outputs/Outputs'
import {
  selectTimeSeriesItemDrawOrderList,
  selectTimeSeriesItemOffset,
  selectTimeSeriesItemShowGrid,
  selectTimeSeriesItemShowLine,
  selectTimeSeriesItemShowTickLabels,
  selectTimeSeriesItemSpan,
  selectTimeSeriesItemXrange,
  selectTimeSeriesItemZeroLine,
  selectVisualizeItemHeight,
  selectVisualizeItemWidth,
  selectTimeSeriesItemKeys,
  selectVisualizeSaveFilename,
  selectVisualizeSaveFormat,
} from 'store/slice/VisualizeItem/VisualizeItemSelectors'
import createColormap from 'colormap'
import { setTimeSeriesItemDrawOrderList } from 'store/slice/VisualizeItem/VisualizeItemSlice'

export const TimeSeriesPlot = React.memo(() => {
  const { itemId, filePath: path } = React.useContext(DisplayDataContext)
  const dispatch = useDispatch()
  const isPending = useSelector(selectTimeSeriesDataIsPending(path))
  const isInitialized = useSelector(selectTimeSeriesDataIsInitialized(path))
  const error = useSelector(selectTimeSeriesDataError(path))
  const isFulfilled = useSelector(selectTimeSeriesDataIsFulfilled(path))

  React.useEffect(() => {
    if (!isInitialized) {
      dispatch(getTimeSeriesInitData({ path, itemId }))
    }
  }, [dispatch, isInitialized, path, itemId])

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
  const drawOrderList = useSelector(selectTimeSeriesItemDrawOrderList(itemId))
  const width = useSelector(selectVisualizeItemWidth(itemId))
  const height = useSelector(selectVisualizeItemHeight(itemId))
  const dataKeys = useSelector(selectTimeSeriesItemKeys(itemId))

  const colorScale = createColormap({
    colormap: 'jet',
    nshades: 100, //maxIndex >= 6 ? maxIndex : 6,
    format: 'hex',
    alpha: 1,
  })

  const data = React.useMemo(() => {
    return Object.fromEntries(
      dataKeys.map((key) => {
        let y = Object.values(timeSeriesData[key])
        const i = Number(key) - 1
        const new_i = Math.floor((i % 10) * 10 + i / 10) % 100
        if (drawOrderList.includes(key) && offset) {
          const activeIdx: number = drawOrderList.findIndex((v) => v === key)
          const mean: number = y.reduce((a, b) => a + b) / y.length
          const std: number =
            span *
            Math.sqrt(y.reduce((a, b) => a + Math.pow(b - mean, 2)) / y.length)
          y = y.map((value) => (value - mean) / (std + 1e-10) + activeIdx)
        }

        return [
          key,
          {
            name: key,
            x: dataXrange,
            y: y,
            visible: drawOrderList.includes(key) ? true : 'legendonly',
            line: { color: colorScale[new_i] },
            error_y: {
              type: 'data',
              array:
                !offset && Object.keys(dataStd).includes(key)
                  ? Object.values(dataStd[key])
                  : null,
              visible: true,
            },
          },
        ]
      }),
    )
  }, [
    timeSeriesData,
    drawOrderList,
    offset,
    span,
    colorScale,
    dataStd,
    dataXrange,
    dataKeys,
  ])

  const annotations = React.useMemo(() => {
    return drawOrderList.map((value) => {
      return {
        x: Number(dataXrange[dataXrange.length - 1]) + dataXrange.length / 10,
        y: data[value].y[dataXrange.length - 1],
        xref: 'x',
        yref: 'y',
        text: `cell: ${value}`,
        arrowhead: 1,
        ax: 0,
        ay: -10,
      }
    })
  }, [data, drawOrderList, dataXrange])

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

  const onLegendClick = (event: LegendClickEvent) => {
    const clickNumber = dataKeys[event.curveNumber]

    const newDrawOrderList = drawOrderList.includes(clickNumber)
      ? drawOrderList.filter((value) => value !== clickNumber)
      : [...drawOrderList, clickNumber]

    dispatch(
      setTimeSeriesItemDrawOrderList({
        itemId,
        drawOrderList: newDrawOrderList,
      }),
    )

    // set DisplayNumbers
    if (!drawOrderList.includes(clickNumber)) {
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
