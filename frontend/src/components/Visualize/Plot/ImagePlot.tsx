import React, { useCallback, useEffect } from 'react'
import PlotlyChart from 'react-plotlyjs-ts'
import { useSelector, useDispatch } from 'react-redux'
import { Button, LinearProgress, TextField, Typography } from '@mui/material'

import { twoDimarrayEqualityFn } from 'utils/EqualityUtils'
import { DisplayDataContext } from '../DataContext'

import {
  selectImageDataError,
  selectImageDataIsInitialized,
  selectImageDataIsPending,
  selectImageDataIsFulfilled,
  selectActiveImageData,
  selectRoiData,
  selectImageDataEndIndex,
} from 'store/slice/DisplayData/DisplayDataSelectors'
import {
  getImageData,
  getRoiData,
  getTimeSeriesDataById,
} from 'store/slice/DisplayData/DisplayDataActions'
import {
  selectImageItemShowticklabels,
  selectImageItemZsmooth,
  selectImageItemShowLine,
  selectImageItemShowGrid,
  selectImageItemShowScale,
  selectImageItemColors,
  selectImageItemActiveIndex,
  selectImageItemStartIndex,
  selectImageItemEndIndex,
  selectRoiItemFilePath,
  selectMultiPlotTimeSeriesItemFilepath,
  selectMultiPlotTimeSeriesItemDisplayNumbers,
  selectRoiItemIndex,
  selectImageItemRoiAlpha,
  selectImageItemDuration,
} from 'store/slice/VisualizeItem/VisualizeItemSelectors'
import {
  incrementImageActiveIndex,
  setImageActiveIndex,
  setImageItemDuration,
  setTimeSeriesItemDisplayNumbers,
} from 'store/slice/VisualizeItem/VisualizeItemSlice'
import { RootState } from 'store/store'
import { Datum, LayoutAxis, PlotData } from 'plotly.js'
import createColormap from 'colormap'
import { getFileName } from 'store/slice/FlowElement/FlowElementUtils'
import Slider from '@mui/material/Slider'
import Box from '@mui/material/Box'

export const ImagePlot = React.memo(() => {
  const { filePath: path, itemId } = React.useContext(DisplayDataContext)

  const startIndex = useSelector(selectImageItemStartIndex(itemId))
  const endIndex = useSelector(selectImageItemEndIndex(itemId))
  const isPending = useSelector(selectImageDataIsPending(path))
  const isInitialized = useSelector(selectImageDataIsInitialized(path))
  const isFulfilled = useSelector(selectImageDataIsFulfilled(path))
  const error = useSelector(selectImageDataError(path))

  const roiFilePath = useSelector(selectRoiItemFilePath(itemId))

  const dispatch = useDispatch()
  React.useEffect(() => {
    if (!isInitialized) {
      dispatch(
        getImageData({
          path,
          startIndex: startIndex ?? 1,
          endIndex: endIndex ?? 10,
        }),
      )
    }
    if (roiFilePath != null) {
      dispatch(getRoiData({ path: roiFilePath }))
    }
  }, [dispatch, isInitialized, path, startIndex, endIndex, roiFilePath])
  if (isPending) {
    return <LinearProgress />
  } else if (error != null) {
    return <Typography color="error">{error}</Typography>
  } else if (isFulfilled) {
    return <ImagePlotImple />
  } else {
    return null
  }
})

const ImagePlotImple = React.memo(() => {
  const { itemId } = React.useContext(DisplayDataContext)
  const activeIndex = useSelector(selectImageItemActiveIndex(itemId))
  return (
    <>
      <ImagePlotChart activeIndex={activeIndex} />
    </>
  )
})

const ImagePlotChart = React.memo<{
  activeIndex: number
}>(({ activeIndex }) => {
  const { filePath: path, itemId } = React.useContext(DisplayDataContext)
  const imageData = useSelector(
    selectActiveImageData(path, activeIndex),
    imageDataEqualtyFn,
  )
  const roiFilePath = useSelector(selectRoiItemFilePath(itemId))
  const roiData = useSelector(
    (state: RootState) =>
      roiFilePath != null ? selectRoiData(roiFilePath)(state) : [],
    imageDataEqualtyFn,
  )

  const maxSize = useSelector(selectImageDataEndIndex(path))
  const showticklabels = useSelector(selectImageItemShowticklabels(itemId))
  const showline = useSelector(selectImageItemShowLine(itemId))
  const zsmooth = useSelector(selectImageItemZsmooth(itemId))
  const showgrid = useSelector(selectImageItemShowGrid(itemId))
  const showscale = useSelector(selectImageItemShowScale(itemId))
  const colorscale = useSelector(selectImageItemColors(itemId))
  const duration = useSelector(selectImageItemDuration(itemId))

  const timeDataMaxIndex = useSelector(selectRoiItemIndex(itemId, roiFilePath))

  const roiAlpha = useSelector(selectImageItemRoiAlpha(itemId))

  const colorscaleRoi = createColormap({
    colormap: 'jet',
    nshades: 100, //timeDataMaxIndex >= 6 ? timeDataMaxIndex : 6,
    format: 'rgba',
    alpha: 1.0,
  })

  const data = React.useMemo(
    () => [
      {
        z: imageData,
        type: 'heatmap',
        name: 'images',
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
        showscale: showscale,
        zsmooth: zsmooth, // ['best', 'fast', false]
        // showlegend: true,
      },
      {
        z: roiData,
        type: 'heatmap',
        name: 'roi',
        hovertemplate: 'cell id: %{z}',
        colorscale: [...Array(timeDataMaxIndex)].map((_, i) => {
          const new_i = Math.floor(((i % 10) * 10 + i / 10) % 100)
          const offset = i / (timeDataMaxIndex - 1)
          const rgba = colorscaleRoi[new_i]
          const hex = rgba2hex(rgba, roiAlpha)
          return [offset, hex]
        }),
        zmin: 1,
        zmax: timeDataMaxIndex,
        hoverongaps: false,
        zsmooth: false,
        showscale: false,
      },
    ],
    [
      imageData,
      roiData,
      zsmooth,
      showscale,
      colorscale,
      colorscaleRoi,
      timeDataMaxIndex,
      roiAlpha,
    ],
  )

  const layout = React.useMemo(
    () => ({
      title: getFileName(path),
      // width: 600,
      // height: 600,
      margin: {
        t: 30, // top
        l: 120, // left
        b: 30, // bottom
      },
      dragmode: 'pan', //'select',
      xaxis: {
        autorange: true,
        showgrid: showgrid,
        showline: showline,
        zeroline: false,
        autotick: true,
        ticks: '',
        showticklabels: showticklabels,
      },
      yaxis: {
        automargin: true,
        autorange: 'reversed',
        showgrid: showgrid,
        showline: showline,
        zeroline: false,
        autotick: true, // todo
        ticks: '',
        showticklabels: showticklabels, // todo
      },
    }),
    [path, showgrid, showline, showticklabels],
  )

  const config = {
    displayModeBar: true,
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
  }, [plotlyHeight, activeIndex])

  const dispatch = useDispatch()
  const timeSeriesFilePath = useSelector(
    selectMultiPlotTimeSeriesItemFilepath(itemId),
  )
  const displayNumbers = useSelector(
    selectMultiPlotTimeSeriesItemDisplayNumbers(itemId),
  )

  const onClick = (event: any) => {
    const points: PlotDatum = event.points[0]
    if (
      timeSeriesFilePath !== null &&
      displayNumbers !== null &&
      points.curveNumber >= 1 &&
      !displayNumbers.includes(points.z - 1)
    ) {
      const newValue = [...displayNumbers, points.z - 1]
      dispatch(
        setTimeSeriesItemDisplayNumbers({
          itemId,
          displayNumbers: newValue,
        }),
      )
      dispatch(
        getTimeSeriesDataById({
          path: timeSeriesFilePath,
          index: points.z - 1,
        }),
      )
    }
  }

  const onSliderChange = (
    event: Event,
    value: number | number[],
    activeThumb: number,
  ) => {
    if (typeof value === 'number') {
      const newIndex = value - 1
      if (newIndex !== activeIndex) {
        dispatch(setImageActiveIndex({ itemId, activeIndex: newIndex }))
      }
    }
  }

  const intervalRef = React.useRef<null | NodeJS.Timeout>(null)

  useEffect(() => {
    if (intervalRef.current !== null) {
      if (activeIndex >= maxSize) {
        clearInterval(intervalRef.current)
        intervalRef.current = null
      }
    }
  }, [activeIndex, maxSize])

  const onPlayClick = useCallback(() => {
    if (activeIndex >= maxSize) {
      dispatch(setImageActiveIndex({ itemId, activeIndex: 0 }))
    }
    if (maxSize > 1 && intervalRef.current === null) {
      intervalRef.current = setInterval(() => {
        dispatch(incrementImageActiveIndex({ itemId }))
      }, duration)
    }
  }, [activeIndex, maxSize, dispatch, duration, itemId])

  const onPauseClick = () => {
    if (intervalRef.current !== null) {
      clearInterval(intervalRef.current)
      intervalRef.current = null
    }
  }

  const onDurationChange = useCallback(
    (event: React.ChangeEvent<HTMLInputElement>) => {
      const newValue =
        event.target.value === '' ? '' : Number(event.target.value)
      if (typeof newValue === 'number') {
        dispatch(setImageItemDuration({ itemId, duration: newValue }))
      }
    },
    [dispatch, itemId],
  )

  return (
    <div ref={ref}>
      <PlotlyChart
        data={data}
        layout={layout}
        config={config}
        onClick={onClick}
      />
      <Box sx={{ width: '50%' }}>
        <Button variant="outlined" onClick={onPlayClick}>
          Play
        </Button>
        <Button variant="outlined" onClick={onPauseClick}>
          Pause
        </Button>
        duration:
        <TextField
          // error={inputError}
          type="number"
          inputProps={{
            step: 100,
            min: 0,
            max: 1000,
          }}
          InputLabelProps={{
            shrink: true,
          }}
          onChange={onDurationChange}
          value={duration}
          // helperText={inputError ? 'index > 0' : undefined}
        />
        msec
        <Typography>Index: {activeIndex + 1}</Typography>
        <Slider
          aria-label="Index"
          defaultValue={1}
          value={activeIndex + 1}
          valueLabelDisplay="auto"
          step={1}
          marks
          min={1}
          max={maxSize === 1 ? 1 : maxSize + 1}
          onChange={onSliderChange}
        />
      </Box>
    </div>
  )
})

function imageDataEqualtyFn(
  a: number[][] | undefined,
  b: number[][] | undefined,
) {
  if (a != null && b != null) {
    return twoDimarrayEqualityFn(a, b)
  } else {
    return a === undefined && b === undefined
  }
}

interface PlotDatum {
  curveNumber: number
  data: PlotData
  customdata: Datum
  pointIndex: number
  pointNumber: number
  x: Datum
  xaxis: LayoutAxis
  y: Datum
  yaxis: LayoutAxis
  z: number
}

function rgba2hex(rgba: [number, number, number, number], alpha: number) {
  const r = rgba[0]
  const g = rgba[1]
  const b = rgba[2]
  const a = alpha

  var outParts = [
    r.toString(16),
    g.toString(16),
    b.toString(16),
    Math.round(a * 255)
      .toString(16)
      .substring(0, 2),
  ]

  // Pad single-digit output values
  outParts.forEach(function (part, i) {
    if (part.length === 1) {
      outParts[i] = '0' + part
    }
  })

  return '#' + outParts.join('')
}
