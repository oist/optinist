import React, {
  useCallback,
  useEffect,
  useState,
  MouseEvent,
  useRef,
} from 'react'
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
import { styled } from '@mui/material/styles'
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
  selectVisualizeSaveFilename,
  selectVisualizeSaveFormat,
  selectImageItemAlpha,
} from 'store/slice/VisualizeItem/VisualizeItemSelectors'
import {
  incrementImageActiveIndex,
  resetAllOrderList,
  setImageActiveIndex,
  setImageItemDuration,
} from 'store/slice/VisualizeItem/VisualizeItemSlice'
import {
  selectingImageArea,
  setImageItemClikedDataId,
} from 'store/slice/VisualizeItem/VisualizeItemActions'
import { addRoiApi, deleteRoiApi, mergeRoiApi } from 'api/outputs/Outputs'

interface PointClick {
  x: number
  y: number
  z: number
}

const WIDTH_CHARTJS = 321
const INIT_WIDTH_ROI = 30

const initSizeDrag = {
  width: INIT_WIDTH_ROI,
  height: INIT_WIDTH_ROI,
  left: Math.floor((WIDTH_CHARTJS - INIT_WIDTH_ROI) / 2),
  top: Math.floor((WIDTH_CHARTJS - INIT_WIDTH_ROI) / 2),
}

enum PositionDrag {
  'LEFT' = 'LEFT',
  'RIGHT' = 'RIGHT',
  'BOTTOM' = 'BOTTOM',
  'TOP' = 'TOP',
}

const sChart = 320

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

  const [isAddRoi, setIsAddRoi] = useState(false)

  const [roiDataState, setRoiDataState] = useState(roiData)

  const [pointClick, setPointClick] = useState<PointClick[]>([])

  const showticklabels = useSelector(selectImageItemShowticklabels(itemId))
  const showline = useSelector(selectImageItemShowLine(itemId))
  const zsmooth = useSelector(selectImageItemZsmooth(itemId))
  const showgrid = useSelector(selectImageItemShowGrid(itemId))
  const showscale = useSelector(selectImageItemShowScale(itemId))
  const colorscale = useSelector(selectImageItemColors(itemId))
  const alpha = useSelector(selectImageItemAlpha(itemId))
  const timeDataMaxIndex = useSelector(selectRoiItemIndex(itemId, roiFilePath))
  const roiAlpha = useSelector(selectImageItemRoiAlpha(itemId))
  const width = useSelector(selectVisualizeItemWidth(itemId))
  const height = useSelector(selectVisualizeItemHeight(itemId))

  const [sizeDrag, setSizeDrag] = useState(initSizeDrag)

  const [startDragAddRoi, setStartDragAddRoi] = useState(false)
  const [positionDrag, setChangeSize] = useState<PositionDrag | undefined>()

  const outputKey: string = useSelector(
    (state: any) => state.visualaizeItem?.items[itemId]?.roiItem?.outputKey,
  )

  const refPageXSize = useRef(0)
  const refPageYSize = useRef(0)

  const colorscaleRoi = createColormap({
    colormap: 'jet',
    nshades: 100, //timeDataMaxIndex >= 6 ? timeDataMaxIndex : 6,
    format: 'rgba',
    alpha: 1.0,
  })

  useEffect(() => {
    setRoiDataState(roiData)
  }, [roiData])

  useEffect(() => {
    onCancel()
    onCancelAdd()
    //eslint-disable-next-line react-hooks/exhaustive-deps
  }, [outputKey])

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
          const rgb = value.rgb
            .replace(/[^0-9,]/g, '')
            .split(',')
            .map((x) => Number(x))
          const hex = rgba2hex(rgb, alpha)
          return [offset, hex]
        }),
        hoverinfo: 'none',
        hoverongaps: false,
        showscale: showscale,
        zsmooth: zsmooth, // ['best', 'fast', false]
      },
      {
        z: roiDataState,
        type: 'heatmap',
        name: 'roi',
        hoverinfo: 'none',
        colorscale: [...Array(timeDataMaxIndex)].map((_, i) => {
          const new_i = Math.floor(((i % 10) * 10 + i / 10) % 100)
          const offset: number = i / (timeDataMaxIndex - 1)
          const rgba = colorscaleRoi[new_i]
          const hex = rgba2hex(rgba, roiAlpha)
          return [offset, hex]
        }),
        zmin: 0,
        zmax: timeDataMaxIndex,
        hoverongaps: false,
        zsmooth: false,
        showscale: false,
      },
    ],
    [
      imageData,
      roiDataState,
      zsmooth,
      showscale,
      colorscale,
      colorscaleRoi,
      timeDataMaxIndex,
      roiAlpha,
      alpha,
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
      height: height - 130,
      margin: {
        t: 30, // top
        l: 100, // left
        b: 20, // bottom
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
    //eslint-disable-next-line react-hooks/exhaustive-deps
    [showgrid, showline, showticklabels, width, height, selectMode, isAddRoi],
  )

  const saveFileName = useSelector(selectVisualizeSaveFilename(itemId))
  const saveFormat = useSelector(selectVisualizeSaveFormat(itemId))

  const config = {
    displayModeBar: true,
    responsive: true,
    toImageButtonOptions: {
      format: saveFormat,
      filename: saveFileName,
      // scale: number;
      // format: 'png' | 'svg' | 'jpeg' | 'webp';
      // height: number;
      // width: number;
    },
  }

  const onClick = (event: any) => {
    const point: PlotDatum = event.points[0]
    if (point.curveNumber >= 1 && outputKey === 'cell_roi') {
      setSelectRoi({
        x: Number(point.x),
        y: Number(point.y),
        z: Number(point.z),
      })
    }
    if (point.curveNumber >= 1 && point.z > 0) {
      dispatch(
        setImageItemClikedDataId({
          itemId,
          clickedDataId: point.z.toString(),
        }),
      )
    }
  }

  const setSelectRoi = (point: PointClick) => {
    if (!point.z) return
    const newPoints = [...pointClick, point]
    const newRoi = roiDataState.map((roi) => {
      return roi.map((element) => {
        if (newPoints.some((p) => p.z === element)) {
          return 0
        }
        return element
      })
    })
    setPointClick([...pointClick, point])
    setRoiDataState(newRoi)
  }

  const onCancel = () => {
    setPointClick([])
    setRoiDataState(roiData)
  }

  const addRoi = () => {
    setIsAddRoi(true)
  }

  const onCancelAdd = () => {
    setIsAddRoi(false)
    setSizeDrag(initSizeDrag)
    setChangeSize(undefined)
  }

  const onMouseDownDragAddRoi = () => {
    setStartDragAddRoi(true)
  }

  const onMouseUpDragAddRoi = () => {
    setStartDragAddRoi(false)
    setChangeSize(undefined)
  }

  const onMouseDownSize = (position: PositionDrag, event: MouseEvent) => {
    setChangeSize(position)
    refPageXSize.current = event.pageX
    refPageYSize.current = event.pageY
  }

  const onMouseMoveAddRoi = (event: MouseEvent<HTMLDivElement>) => {
    const { pageX, pageY } = event
    let newSizeDrag
    if (startDragAddRoi) {
      const { y } = (event.currentTarget as any).getBoundingClientRect()
      let newX = sizeDrag.left + (pageX - refPageXSize.current)
      let newY = Math.ceil(pageY - y - 15)
      if (newX < 0) newX = 0
      else if (newX + sizeDrag.width > sChart) newX = sChart - sizeDrag.width
      if (newY < 0) newY = 0
      else if (newY + sizeDrag.height > sChart) newY = sChart - sizeDrag.height
      newSizeDrag = { ...sizeDrag, left: newX, top: newY }
    } else if (positionDrag === PositionDrag.LEFT) {
      const newWidth = sizeDrag.width - (pageX - refPageXSize.current)
      const newLeft = sizeDrag.left + (pageX - refPageXSize.current)
      if (newWidth < 10 || newLeft < 1) return
      newSizeDrag = { ...sizeDrag, width: newWidth, left: newLeft }
    } else if (positionDrag === PositionDrag.RIGHT) {
      const newWidth = sizeDrag.width + (pageX - refPageXSize.current)
      if (newWidth < 10 || newWidth > sChart - sizeDrag.left) return
      newSizeDrag = { ...sizeDrag, width: newWidth }
    } else if (positionDrag === PositionDrag.BOTTOM) {
      const newHeight = sizeDrag.height + (pageY - refPageYSize.current)
      if (newHeight < 10 || newHeight > sChart - sizeDrag.top) return
      newSizeDrag = { ...sizeDrag, height: newHeight }
    } else if (positionDrag === PositionDrag.TOP) {
      const newHeight = sizeDrag.height - (pageY - refPageYSize.current)
      const newTop = sizeDrag.top + (pageY - refPageYSize.current)
      if (newHeight < 10 || newTop < 1) return
      newSizeDrag = { ...sizeDrag, height: newHeight, top: newTop }
    }
    if (newSizeDrag) setSizeDrag({ ...sizeDrag, ...newSizeDrag })
    refPageXSize.current = pageX
    refPageYSize.current = pageY
  }

  const addRoiSubmit = async () => {
    if (!roiFilePath) return
    const sizeX = roiDataState[0].length - 1
    const sizeY = roiDataState.length - 1
    const xAdd = Number(((sizeDrag.width + 2) / (sChart / sizeX)).toFixed(1))
    const yAdd = Number(((sizeDrag.height + 2) / (sChart / sizeY)).toFixed(1))
    const x = Number((sizeDrag.left / (sChart / sizeX)).toFixed(1))
    const y = Number((sizeDrag.top / (sChart / sizeY)).toFixed(1))

    const pointCenter = {
      posx: x + Math.floor(xAdd / 2),
      posy: y + Math.floor(yAdd / 2),
      sizex: xAdd,
      sizey: yAdd,
    }
    try {
      await addRoiApi(roiFilePath, pointCenter)
    } catch {}
    onCancelAdd()
    dispatch(getRoiData({ path: roiFilePath }))
  }

  const onMergeRoi = async () => {
    if (!roiFilePath) return
    dispatch(resetAllOrderList())
    try {
      await mergeRoiApi(roiFilePath, {
        ids: pointClick.map((point) => point.z - 1),
      })
    } catch {}
    onCancel()
    dispatch(getRoiData({ path: roiFilePath }))
  }

  const onDeleteRoi = async () => {
    if (!roiFilePath) return
    try {
      await deleteRoiApi(roiFilePath, {
        ids: pointClick.map((point) => point.z - 1),
      })
    } catch {}
    onCancel()
    dispatch(getRoiData({ path: roiFilePath }))
  }

  const renderActionRoi = () => {
    if (!roiDataState?.length || outputKey !== 'cell_roi') return null
    if (!isAddRoi) {
      return <LinkDiv onClick={addRoi}>Add Roi</LinkDiv>
    }
    return (
      <BoxDiv>
        <LinkDiv onClick={addRoiSubmit}>Ok</LinkDiv>
        <LinkDiv onClick={onCancelAdd}>Cancel</LinkDiv>
      </BoxDiv>
    )
  }

  return (
    <div>
      <Box sx={{ display: 'flex' }}>
        <Box sx={{ flexGrow: 1, mt: 1 }}>
          <PlayBack activeIndex={activeIndex} />
        </Box>
        <FormControlLabel
          sx={{ ml: 1 }}
          control={<Switch checked={selectMode} onChange={handleChange} />}
          label="drag select"
        />
      </Box>
      <Box sx={{ minHeight: 5.5 }}>
        {pointClick.length ? (
          <>
            <BoxDiv>
              <span>Roi Selecteds: [{String(pointClick.map((e) => e.z))}]</span>
            </BoxDiv>
            <BoxDiv>
              {pointClick.length >= 2 ? (
                <LinkDiv sx={{ ml: 0 }} onClick={onMergeRoi}>
                  Merge Roi
                </LinkDiv>
              ) : null}
              <LinkDiv sx={{ color: '#F84E1B' }} onClick={onDeleteRoi}>
                Delete Roi
              </LinkDiv>
              <LinkDiv onClick={onCancel}>Cancel</LinkDiv>
            </BoxDiv>
          </>
        ) : (
          renderActionRoi()
        )}
      </Box>
      <div style={{ position: 'relative' }}>
        <PlotlyChart
          data={data}
          layout={layout}
          config={config}
          onClick={onClick}
          onSelecting={onSelecting}
        />
        {isAddRoi ? (
          <DivAddRoi>
            <DivSvg
              onMouseLeave={onMouseUpDragAddRoi}
              onMouseMove={onMouseMoveAddRoi}
              onMouseUp={onMouseUpDragAddRoi}
            >
              <DivDrag style={sizeDrag}>
                <DragCenter
                  onMouseDown={onMouseDownDragAddRoi}
                  style={{
                    width: sizeDrag.width - 1,
                    height: sizeDrag.height - 1,
                    cursor: !startDragAddRoi ? 'grab' : 'grabbing',
                  }}
                />
                <DragSizeLeft
                  onMouseDown={(event) =>
                    onMouseDownSize(PositionDrag.LEFT, event)
                  }
                />
                <DragSizeRight
                  onMouseDown={(event) => {
                    onMouseDownSize(PositionDrag.RIGHT, event)
                  }}
                />
                <DragSizeTop
                  onMouseDown={(event) => {
                    onMouseDownSize(PositionDrag.TOP, event)
                  }}
                />
                <DragSizeBottom
                  onMouseDown={(event) => {
                    onMouseDownSize(PositionDrag.BOTTOM, event)
                  }}
                />
              </DivDrag>
            </DivSvg>
          </DivAddRoi>
        ) : null}
      </div>
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
      <Button sx={{ mt: 1.5 }} variant="outlined" onClick={onPlayClick}>
        Play
      </Button>
      <Button sx={{ mt: 1.5, ml: 1 }} variant="outlined" onClick={onPauseClick}>
        Pause
      </Button>
      <TextField
        sx={{ width: 100, ml: 2 }}
        label="Duration [msec]"
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

function rgba2hex(rgba: number[], alpha: number) {
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

  return `#${outParts.join('')}`
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

const BoxDiv = styled('div')({
  mt: 1,
  display: 'flex',
  alignItems: 'center',
  listStyle: 'none',
  padding: 0,
  margin: 0,
})

const LinkDiv = styled('div')({
  marginLeft: 16,
  textDecoration: 'underline',
  cursor: 'pointer',
  color: '#1155cc',
  zIndex: 999,
  position: 'relative',
})

const DivAddRoi = styled('div')({
  width: '100%',
  height: '100%',
  position: 'absolute',
  left: 0,
  top: 0,
  borderRadius: 100,
})

const DivSvg = styled('div')({
  width: 321,
  height: 321,
  marginTop: 30,
  marginLeft: 99,
  position: 'relative',
})

const DivDrag = styled('div')({
  border: '1px solid #ffffff',
  position: 'absolute',
  borderRadius: 100,
})

const DragCenter = styled('div')({
  borderRadius: 100,
  cursor: 'grab',
})

const DragSize = styled('div')({
  width: 3,
  height: 3,
  borderRadius: 100,
  position: 'absolute',
  background: '#fff',
})

const DragSizeLeft = styled(DragSize)({
  top: `calc(50% - 1px)`,
  left: -2,
  cursor: 'ew-resize',
})

const DragSizeRight = styled(DragSize)({
  top: `calc(50% - 1px)`,
  right: -2,
  cursor: 'ew-resize',
})

const DragSizeTop = styled(DragSize)({
  top: -2,
  right: `calc(50% - 1px)`,
  cursor: 'ns-resize',
})

const DragSizeBottom = styled(DragSize)({
  bottom: -2,
  right: `calc(50% - 1px)`,
  cursor: 'ns-resize',
})
