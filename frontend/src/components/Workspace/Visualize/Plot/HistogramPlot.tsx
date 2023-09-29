import React from 'react'
import { useSelector, useDispatch } from 'react-redux'
import PlotlyChart from 'react-plotlyjs-ts'
import { DisplayDataContext } from '../DataContext'
import {
  selectHistogramData,
  selectHistogramDataError,
  selectHistogramDataIsFulfilled,
  selectHistogramDataIsInitialized,
  selectHistogramDataIsPending,
  selectHistogramMeta,
} from 'store/slice/DisplayData/DisplayDataSelectors'
import { getHistogramData } from 'store/slice/DisplayData/DisplayDataActions'
import {
  Box,
  FormControl,
  Input,
  InputLabel,
  LinearProgress,
  Typography,
} from '@mui/material'
import {
  selectHistogramItemBins,
  selectVisualizeItemHeight,
  selectVisualizeItemWidth,
} from 'store/slice/VisualizeItem/VisualizeItemSelectors'
import { setHistogramItemBins } from 'store/slice/VisualizeItem/VisualizeItemSlice'

export const HistogramPlot = React.memo(() => {
  const { filePath: path } = React.useContext(DisplayDataContext)
  const dispatch = useDispatch()
  const isPending = useSelector(selectHistogramDataIsPending(path))
  const isInitialized = useSelector(selectHistogramDataIsInitialized(path))
  const error = useSelector(selectHistogramDataError(path))
  const isFulfilled = useSelector(selectHistogramDataIsFulfilled(path))

  React.useEffect(() => {
    if (!isInitialized) {
      dispatch(getHistogramData({ path }))
    }
  }, [dispatch, isInitialized, path])

  if (isPending) {
    return <LinearProgress />
  } else if (error != null) {
    return <Typography color="error">{error}</Typography>
  } else if (isFulfilled) {
    return <HistogramPlotImple />
  } else {
    return null
  }
})

const HistogramPlotImple = React.memo(() => {
  const { filePath: path, itemId } = React.useContext(DisplayDataContext)
  const histogramData = useSelector(selectHistogramData(path))
  const meta = useSelector(selectHistogramMeta(path))
  const width = useSelector(selectVisualizeItemWidth(itemId))
  const height = useSelector(selectVisualizeItemHeight(itemId))
  const bins = useSelector(selectHistogramItemBins(itemId))

  const data = React.useMemo(
    () =>
      histogramData != null
        ? [
            {
              x: histogramData[0],
              type: 'histogram',
              autobinx: false,
              xbins: {
                size:
                  (Math.max(...histogramData[0]) -
                    Math.min(...histogramData[0])) /
                  (bins - 1),
              },
            },
          ]
        : [],
    [histogramData, bins],
  )

  const layout = React.useMemo(
    () => ({
      title: {
        text: meta?.title,
        x: 0.1,
      },
      width: width,
      height: height - 120,
      dragmode: 'pan',
      margin: {
        t: 50, // top
        l: 50, // left
        b: 40, // bottom
      },
      autosize: true,
      xaxis: {
        title: meta?.xlabel,
      },
      yaxis: {
        title: meta?.ylabel,
      },
    }),
    [meta, width, height],
  )

  return (
    <div>
      <Box sx={{ display: 'flex' }}>
        <Box sx={{ flexGrow: 1, ml: 1 }}>
          <InputBins />
        </Box>
      </Box>
      <PlotlyChart data={data} layout={layout} />
    </div>
  )
})

const InputBins = React.memo(() => {
  const { itemId } = React.useContext(DisplayDataContext)
  const dispatch = useDispatch()
  const bins = useSelector(selectHistogramItemBins(itemId))

  const handleChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    dispatch(
      setHistogramItemBins({ itemId, bins: parseInt(event.target.value) }),
    )
  }

  return (
    <FormControl variant="standard" sx={{ m: 1, minWidth: 120 }}>
      <InputLabel>bins</InputLabel>
      <Input type="number" value={bins} onChange={handleChange} />
    </FormControl>
  )
})
