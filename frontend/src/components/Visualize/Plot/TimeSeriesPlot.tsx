import React, { useEffect } from 'react'
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
  selectTimeSeriesItemShowGrid,
  selectTimeSeriesItemShowLine,
  selectTimeSeriesItemShowTickLabels,
  selectTimeSeriesItemSpan,
  selectTimeSeriesItemXrange,
  selectTimeSeriesItemZeroLine,
} from 'store/slice/VisualizeItem/VisualizeItemSelectors'
import { setTimeSeriesItemCheckedList } from 'store/slice/VisualizeItem/VisualizeItemSlice'
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
  const displayNumbers = useSelector(selectTimeSeriesItemDisplayNumbers(itemId))
  const checkedList = useSelector(selectTimeSeriesItemCheckedList(itemId))

  const colorScale = createColormap({
    colormap: 'jet',
    nshades: 100, //maxIndex >= 6 ? maxIndex : 6,
    format: 'hex',
    alpha: 1,
  })

  React.useEffect(() => {
    const keys = Object.keys(timeSeriesData)
    if (checkedList.length === 0 && keys.length !== 0) {
      const checkedList = keys.map((_, i) => {
        if (i === 0) {
          return true
        }
        return false
      })
      dispatch(
        setTimeSeriesItemCheckedList({
          itemId,
          checkedList,
        }),
      )
    }
  }, [timeSeriesData])

  const data = React.useMemo(() => {
    if (timeSeriesData === null) {
      return []
    }
    return Object.keys(timeSeriesData).map((key, i) => {
      let y = Object.values(timeSeriesData[key])
      const new_i = Math.floor((i % 10) * 10 + i / 10) % 100

      if (displayNumbers.includes(i)) {
        if (offset) {
          const activeIdx: number = displayNumbers.findIndex(
            (value) => value === i,
          )
          const mean: number = y.reduce((a, b) => a + b) / y.length
          const std: number =
            span *
            Math.sqrt(y.reduce((a, b) => a + Math.pow(b - mean, 2)) / y.length)
          return {
            name: `(${String(parseInt(key) + 1)})`,
            x: dataXrange,
            y: y.map((value) => (value - mean) / std + activeIdx),
            visible: true,
            line: { color: colorScale[new_i] },
            error_y: {
              type: 'data',
              array: null,
              visible: true,
            },
          }
        } else {
          return {
            name: `(${String(parseInt(key) + 1)})`,
            x: dataXrange,
            y: y,
            visible: true,
            line: { color: colorScale[new_i] },
            error_y: {
              type: 'data',
              array: Object.keys(dataStd).includes(key)
                ? Object.values(dataStd[key])
                : null,
              visible: true,
            },
          }
        }
      } else {
        return {
          name: `(${String(parseInt(key) + 1)})`,
          x: dataXrange,
          y: y,
          visible: 'legendonly',
          line: { color: colorScale[new_i] },
          error_y: {
            type: 'data',
            array: Object.keys(dataStd).includes(key)
              ? Object.values(dataStd[key])
              : null,
            visible: true,
          },
        }
      }
    })
  }, [timeSeriesData, displayNumbers, offset, span, colorScale])

  const getAnnotation = () => {
    if (data.length !== 0) {
      return displayNumbers.map((i) => {
        return {
          x: dataXrange[0],
          y: Math.max(...data[i].y),
          xref: 'x',
          yref: 'y',
          text: `cell: ${i + 1}`,
          arrowhead: 1,
          ax: 0,
          ay: -10,
        }
      })
    } else {
      return []
    }
  }

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
      annotations: getAnnotation(),
    }),
    [
      path,
      xrange,
      showgrid,
      showline,
      showticklabels,
      zeroline,
      offset,
      span,
      data,
    ],
  )

  const config = {
    displayModeBar: true,
    responsive: true,
  }

  const onLegendClick = (event: LegendClickEvent) => {
    const clickNumber = event.curveNumber

    // set DisplayNumbers
    if (displayNumbers.includes(clickNumber)) {
      dispatch(
        setTimeSeriesItemDisplayNumbers({
          itemId,
          displayNumbers: displayNumbers.filter(
            (value) => value !== clickNumber,
          ),
        }),
      )

      dispatch(
        setTimeSeriesItemCheckedList({
          itemId,
          checkedList: checkedList.map((v, i) => {
            if (i === clickNumber) {
              return false
            }
            return v
          }),
        }),
      )
    } else {
      dispatch(
        setTimeSeriesItemDisplayNumbers({
          itemId,
          displayNumbers: [...displayNumbers, clickNumber],
        }),
      )

      dispatch(
        setTimeSeriesItemCheckedList({
          itemId,
          checkedList: checkedList.map((v, i) => {
            if (i === clickNumber) {
              return true
            }
            return v
          }),
        }),
      )

      dispatch(getTimeSeriesDataById({ path, index: clickNumber }))
    }

    return false
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
      <PlotlyChart
        data={data}
        layout={layout}
        config={config}
        onLegendClick={onLegendClick}
      />
    </div>
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
