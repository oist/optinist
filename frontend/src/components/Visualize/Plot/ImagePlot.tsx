import React, { useCallback, useEffect } from 'react'
import PlotlyChart from 'react-plotlyjs-ts'
import { useSelector, useDispatch } from 'react-redux'
import { RootState } from 'store/store'
import { Datum, LayoutAxis, PlotData, PlotSelectionEvent } from 'plotly.js'
import createColormap from 'colormap'
import { Button, LinearProgress, TextField, Typography } from '@mui/material'
import FormControlLabel from '@mui/material/FormControlLabel'
import Switch from '@mui/material/Switch'
import Slider from '@mui/material/Slider'
import Box from '@mui/material/Box'

import { twoDimarrayEqualityFn } from 'utils/EqualityUtils'
import { DisplayDataContext } from '../DataContext'

import {
  selectImageDataError,
  selectImageDataIsInitialized,
  selectImageDataIsPending,
  selectImageDataIsFulfilled,
  selectActiveImageData,
  selectRoiData,
  selectImageDataMaxSize,
} from 'store/slice/DisplayData/DisplayDataSelectors'
import {
  getImageData,
  getRoiData,
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
  selectRoiItemIndex,
  selectImageItemRoiAlpha,
  selectImageItemDuration,
  selectVisualizeItemWidth,
  selectVisualizeItemHeight,
} from 'store/slice/VisualizeItem/VisualizeItemSelectors'
import {
  incrementImageActiveIndex,
  setImageActiveIndex,
  setImageItemDuration,
} from 'store/slice/VisualizeItem/VisualizeItemSlice'
import {
  selectingImageArea,
  setImageItemClikedDataId,
} from 'store/slice/VisualizeItem/VisualizeItemActions'

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
  return <ImagePlotChart activeIndex={activeIndex} />
})

const ImagePlotChart = React.memo<{
  activeIndex: number
}>(({ activeIndex }) => {
  const dispatch = useDispatch()
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

  const showticklabels = useSelector(selectImageItemShowticklabels(itemId))
  const showline = useSelector(selectImageItemShowLine(itemId))
  const zsmooth = useSelector(selectImageItemZsmooth(itemId))
  const showgrid = useSelector(selectImageItemShowGrid(itemId))
  const showscale = useSelector(selectImageItemShowScale(itemId))
  const colorscale = useSelector(selectImageItemColors(itemId))
  const timeDataMaxIndex = useSelector(selectRoiItemIndex(itemId, roiFilePath))
  const roiAlpha = useSelector(selectImageItemRoiAlpha(itemId))
  const width = useSelector(selectVisualizeItemWidth(itemId))
  const height = useSelector(selectVisualizeItemHeight(itemId))

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

  const [selectMode, setSelectMode] = React.useState(false)

  const handleChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    setSelectMode(event.target.checked)
  }
  // debounceでイベントを間引きする。onSelectedはそれっぽい名前だが動かなかった。
  const onSelecting = debounce((event: PlotSelectionEvent) => {
    if (event.range != null) {
      dispatch(selectingImageArea({ itemId, range: event.range }))
    }
  })
  const layout = React.useMemo(
    () => ({
      width: width,
      height: height - 120 - 9000 / width,
      margin: {
        t: 30, // top
        l: 120, // left
        b: 30, // bottom
      },
      dragmode: selectMode ? 'select' : 'pan',
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
    [showgrid, showline, showticklabels, width, height, selectMode],
  )

  const config = {
    displayModeBar: true,
    responsive: true,
  }

  const onClick = (event: any) => {
    const points: PlotDatum = event.points[0]
    if (points.curveNumber >= 1) {
      dispatch(
        setImageItemClikedDataId({ itemId, clickedDataId: points.z - 1 }),
      )
    }
  }

  return (
    <div>
      <Box sx={{ display: 'flex' }}>
        <Box sx={{ flexGrow: 1, ml: 1 }}>
          <PlayBack activeIndex={activeIndex} />
        </Box>
        <FormControlLabel
          sx={{ m: 1 }}
          control={<Switch checked={selectMode} onChange={handleChange} />}
          label="drag select"
        />
      </Box>
      <PlotlyChart
        data={data}
        layout={layout}
        config={config}
        onClick={onClick}
        onSelecting={onSelecting}
      />
    </div>
  )
})

const PlayBack = React.memo<{ activeIndex: number }>(({ activeIndex }) => {
  const dispatch = useDispatch()
  const { filePath: path, itemId } = React.useContext(DisplayDataContext)

  const maxSize = useSelector(selectImageDataMaxSize(path))
  const startIndex = useSelector(selectImageItemStartIndex(itemId))
  const endIndex = useSelector(selectImageItemEndIndex(itemId))
  const duration = useSelector(selectImageItemDuration(itemId))

  const onSliderChange = (
    event: Event,
    value: number | number[],
    activeThumb: number,
  ) => {
    if (typeof value === 'number') {
      const newIndex = value - startIndex
      if (newIndex >= 0 && newIndex !== activeIndex) {
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
    <>
      <Button variant="outlined" onClick={onPlayClick}>
        Play
      </Button>
      <Button variant="outlined" onClick={onPauseClick}>
        Pause
      </Button>
      <TextField
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
      />
      <Slider
        aria-label="Custom marks"
        defaultValue={20}
        value={startIndex + activeIndex}
        valueLabelDisplay="auto"
        step={1}
        marks
        min={startIndex}
        max={maxSize === 0 ? 0 : endIndex}
        onChange={onSliderChange}
      />
    </>
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

function debounce<T extends (...args: any[]) => unknown>(
  callback: T,
  delay = 500,
): (...args: Parameters<T>) => void {
  let timeoutId: NodeJS.Timeout
  return (...args) => {
    clearTimeout(timeoutId)
    timeoutId = setTimeout(() => callback(...args), delay)
  }
}
