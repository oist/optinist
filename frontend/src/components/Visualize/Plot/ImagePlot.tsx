import React from 'react'
import PlotlyChart from 'react-plotlyjs-ts'
import { useSelector, useDispatch } from 'react-redux'
import {
  Button,
  LinearProgress,
  MobileStepper,
  Typography,
  useTheme,
} from '@material-ui/core'
import KeyboardArrowLeft from '@material-ui/icons/KeyboardArrowLeft'
import KeyboardArrowRight from '@material-ui/icons/KeyboardArrowRight'

import { twoDimarrayEqualityFn } from 'utils/EqualityUtils'
import { DisplayDataContext } from '../DisplayDataItem'

import {
  selectImageDataError,
  selectImageDataIsInitialized,
  selectImageDataIsPending,
  selectImageDataMaxIndex,
  selectImageDataIsFulfilled,
  selectActiveImageData,
} from 'store/slice/DisplayData/DisplayDataSelectors'
import { getImageData } from 'store/slice/DisplayData/DisplayDataActions'
import { selectImageMaxIndexByNodeId } from 'store/slice/InputNode/InputNodeSelectors'
import {
  selectImageItemShowticklabels,
  selectImageItemZsmooth,
  selectImageItemShowLine,
  selectImageItemShowGrid,
  selectImageItemShowScale,
  selectImageItemColors,
  selectImageItemActiveIndex,
} from 'store/slice/VisualizeItem/VisualizeItemSelectors'
import { RootState } from 'store/store'
import {
  decrementImageActiveIndex,
  incrementImageActiveIndex,
} from 'store/slice/VisualizeItem/VisualizeItemSlice'

export const ImagePlot = React.memo(() => {
  const { filePath: path, nodeId } = React.useContext(DisplayDataContext)
  const maxIndex = useSelector((state: RootState) => {
    if (nodeId) {
      return selectImageMaxIndexByNodeId(nodeId)(state)
    } else {
      return 1
    }
  })
  const dispatch = useDispatch()
  const isPending = useSelector(selectImageDataIsPending(path))
  const isInitialized = useSelector(selectImageDataIsInitialized(path))
  const isFulfilled = useSelector(selectImageDataIsFulfilled(path))
  const error = useSelector(selectImageDataError(path))
  React.useEffect(() => {
    if (!isInitialized) {
      dispatch(getImageData({ path, maxIndex: maxIndex ?? 1 }))
    }
  }, [dispatch, isInitialized, path, maxIndex])
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
  const { filePath: path, itemId } = React.useContext(DisplayDataContext)
  const maxIndex = useSelector(selectImageDataMaxIndex(path))
  const activeIndex = useSelector(selectImageItemActiveIndex(itemId))
  const dispatch = useDispatch()
  const handleNext = () =>
    dispatch(incrementImageActiveIndex({ itemId, activeIndex }))
  const handleBack = () =>
    dispatch(decrementImageActiveIndex({ itemId, activeIndex }))
  const theme = useTheme()
  return (
    <>
      <MobileStepper
        steps={(maxIndex ?? 0) + 1}
        position="static"
        variant="text"
        activeStep={activeIndex}
        nextButton={
          <Button
            size="small"
            onClick={handleNext}
            disabled={activeIndex === (maxIndex ?? 0)}
          >
            <Typography>Next</Typography>
            {theme.direction === 'rtl' ? (
              <KeyboardArrowLeft />
            ) : (
              <KeyboardArrowRight />
            )}
          </Button>
        }
        backButton={
          <Button
            size="small"
            onClick={handleBack}
            disabled={activeIndex === 0}
          >
            {theme.direction === 'rtl' ? (
              <KeyboardArrowRight />
            ) : (
              <KeyboardArrowLeft />
            )}
            <Typography>Back</Typography>
          </Button>
        }
      />
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
  // const testData1 = [
  //   [0, 10, 30],
  //   [10, 20, 10],
  //   [0, 10, 10],
  //   [20, 0, 10],
  // ]
  // const testData2 = [
  //   [null, 1, 2],
  //   [3, null, 5],
  //   [6, 7, 8],
  //   [9, 10, 11],
  // ]
  // const data = React.useMemo(
  //   () => [
  //     {
  //       z: testData1,
  //       type: 'heatmap',
  //       name: 'images',
  //       colorscale: [
  //         [0, '#000000'],
  //         [1, '#ffffff'],
  //       ],
  //       hoverongaps: false,
  //       showscale: false,
  //       zsmooth: false,   // ['best', 'fast', false]
  //       showlegend: true,
  //     },
  //     {
  //       z: testData2,
  //       type: 'heatmap',
  //       name: 'iscell',
  //       colorscale:  [
  //         [0, 'rgb(166,206,227)'],
  //         [0.25, 'rgb(31,120,180)'],
  //         [0.45, 'rgb(178,223,138)'],
  //         [0.65, 'rgb(51,160,44)'],
  //         [0.85, 'rgb(251,154,153)'],
  //         [1, 'rgb(227,26,28)'],
  //       ],
  //       hoverongaps: false,
  //       showscale: false,
  //       showlegend: true
  //     },
  //   ],
  //   [testData1, testData2],
  // )

  const showticklabels = useSelector(selectImageItemShowticklabels(itemId))
  const showline = useSelector(selectImageItemShowLine(itemId))
  const zsmooth = useSelector(selectImageItemZsmooth(itemId))
  const showgrid = useSelector(selectImageItemShowGrid(itemId))
  const showscale = useSelector(selectImageItemShowScale(itemId))
  const colorscale = useSelector(selectImageItemColors(itemId))

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
        showlegend: true,
      },
    ],
    [imageData, zsmooth, showscale, colorscale],
  )

  const layout = {
    title: path.split('/').reverse()[0],
    margin: {
      t: 30, // top
      l: 120, // left
      b: 30, // bottom
    },
    dragmode: 'pan',
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
      autorange: 'reversed',
      showgrid: showgrid,
      showline: showline,
      zeroline: false,
      autotick: true, // todo
      ticks: '',
      showticklabels: showticklabels, // todo
    },
  }
  const config = {
    displayModeBar: true,
    scrollZoom: true,
  }
  return (
    <div className="imagePlotChart">
      <PlotlyChart data={data} layout={layout} config={config} />
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
