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
  decrementImageActiveIndex,
  incrementImageActiveIndex,
} from 'store/slice/DisplayData/DisplayDataSlice'
import {
  selectActiveImageData,
  selectImageDataActiveIndex,
  selectImageDataError,
  selectImageDataIsInitialized,
  selectImageDataIsPending,
  selectImageDataMaxIndex,
  selectImageDataIsFulfilled,
} from 'store/slice/DisplayData/DisplayDataSelectors'
import { getImageData } from 'store/slice/DisplayData/DisplayDataActions'
import { selectImageMaxIndexByNodeId } from 'store/slice/InputNode/InputNodeSelectors'
import { selectNodeLabelById } from 'store/slice/FlowElement/FlowElementSelectors'

export const ImagePlot = React.memo(() => {
  const { filePath: path, nodeId } = React.useContext(DisplayDataContext)
  const maxIndex = useSelector(selectImageMaxIndexByNodeId(nodeId))
  const dispatch = useDispatch()
  const isPending = useSelector(selectImageDataIsPending(path))
  const isInitialized = useSelector(selectImageDataIsInitialized(path))
  const isFulfilled = useSelector(selectImageDataIsFulfilled(path))
  const error = useSelector(selectImageDataError(path))
  React.useEffect(() => {
    if (!isInitialized) {
      dispatch(getImageData({ path, maxIndex: maxIndex ?? undefined }))
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
  const { filePath: path } = React.useContext(DisplayDataContext)
  const maxIndex = useSelector(selectImageDataMaxIndex(path))
  const activeIndex = useSelector(selectImageDataActiveIndex(path))
  const dispatch = useDispatch()
  const handleNext = () => dispatch(incrementImageActiveIndex({ path }))
  const handleBack = () => dispatch(decrementImageActiveIndex({ path }))
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
      <ImagePlotChart />
    </>
  )
})

const ImagePlotChart = React.memo(() => {
  const { filePath: path, nodeId } = React.useContext(DisplayDataContext)
  const label = useSelector(selectNodeLabelById(nodeId))
  const imageData = useSelector(selectActiveImageData(path), imageDataEqualtyFn)
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

  const data = React.useMemo(
    () => [
      {
        z: imageData,
        type: 'heatmap',
        name: 'images',
        colorscale: [
          // todo グラデーションの設定
          [0, '#000000'],
          [1, '#ffffff'],
        ],
        hoverongaps: false,
        showscale: false,
        zsmooth: 'best', // ['best', 'fast', false]
        showlegend: true,
      },
    ],
    [imageData],
  )

  const layout = {
    title: label,
    margin: {
      t: 30, // top
      l: 120, // left
      b: 30, // bottom
    },
    dragmode: 'pan',
    xaxis: {
      autorange: true,
      showgrid: false,
      zeroline: false,
      showline: false,
      autotick: true,
      ticks: '',
      showticklabels: false,
    },
    yaxis: {
      autorange: 'reversed',
      showgrid: false,
      zeroline: false,
      showline: false,
      autotick: true, // todo
      ticks: '',
      showticklabels: true, // todo
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
