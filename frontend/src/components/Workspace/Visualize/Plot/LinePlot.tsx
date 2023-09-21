import React from 'react'
import { useSelector, useDispatch } from 'react-redux'
import PlotlyChart from 'react-plotlyjs-ts'
import { DisplayDataContext } from '../DataContext'
import {
  selectLineColumns,
  selectLineData,
  selectLineDataError,
  selectLineDataIsFulfilled,
  selectLineDataIsInitialized,
  selectLineDataIsPending,
  selectLineIndex,
} from 'store/slice/DisplayData/DisplayDataSelectors'
import { getLineData } from 'store/slice/DisplayData/DisplayDataActions'
import {
  Box,
  FormControl,
  InputLabel,
  LinearProgress,
  MenuItem,
  Select,
  SelectChangeEvent,
  Typography,
} from '@mui/material'
import {
  selectLineItemSelectedIndex,
  selectVisualizeItemHeight,
  selectVisualizeItemWidth,
} from 'store/slice/VisualizeItem/VisualizeItemSelectors'
import { setLineItemSelectedIndex } from 'store/slice/VisualizeItem/VisualizeItemSlice'
export const LinePlot = React.memo(() => {
  const { filePath: path } = React.useContext(DisplayDataContext)
  const dispatch = useDispatch()
  const isPending = useSelector(selectLineDataIsPending(path))
  const isInitialized = useSelector(selectLineDataIsInitialized(path))
  const error = useSelector(selectLineDataError(path))
  const isFulfilled = useSelector(selectLineDataIsFulfilled(path))

  React.useEffect(() => {
    if (!isInitialized) {
      dispatch(getLineData({ path }))
    }
  }, [dispatch, isInitialized, path])

  if (isPending) {
    return <LinearProgress />
  } else if (error != null) {
    return <Typography color="error">{error}</Typography>
  } else if (isFulfilled) {
    return <LinePlotImple />
  } else {
    return null
  }
})

const LinePlotImple = React.memo(() => {
  const { filePath: path, itemId } = React.useContext(DisplayDataContext)
  const lineData = useSelector(selectLineData(path))
  const columns = useSelector(selectLineColumns(path))
  const index = useSelector(selectLineIndex(path))
  const selectedIndex = useSelector(selectLineItemSelectedIndex(itemId))
  const width = useSelector(selectVisualizeItemWidth(itemId))
  const height = useSelector(selectVisualizeItemHeight(itemId))

  const data = React.useMemo(
    () =>
      lineData != null
        ? [
            {
              x: columns,
              y: lineData[selectedIndex],
            },
          ]
        : [],
    [lineData, columns, selectedIndex],
  )

  const layout = React.useMemo(
    () => ({
      width: width,
      height: height - 120,
      dragmode: 'pan',
      margin: {
        t: 60, // top
        l: 50, // left
        b: 30, // bottom
      },
      autosize: true,
      xaxis: {
        tickvals: columns,
      },
    }),
    [width, height, columns],
  )

  return (
    <>
      <Box sx={{ display: 'flex' }}>
        <Box sx={{ flexGrow: 1, ml: 1 }}>
          <LineItemIndexSelect index={index} />
        </Box>
      </Box>
      <PlotlyChart data={data} layout={layout} />
    </>
  )
})

const LineItemIndexSelect = React.memo<{ index: number[] }>(({ index }) => {
  const { itemId } = React.useContext(DisplayDataContext)
  const dispatch = useDispatch()
  const selectedIndex = useSelector(selectLineItemSelectedIndex(itemId))

  const handleChange = (event: SelectChangeEvent<number>) => {
    dispatch(
      setLineItemSelectedIndex({
        itemId,
        selectedIndex: Number(event.target.value),
      }),
    )
  }

  return (
    <FormControl variant="standard" sx={{ m: 1, minWidth: 120 }}>
      <InputLabel>index</InputLabel>
      <Select label="smooth" value={selectedIndex} onChange={handleChange}>
        {index.map((_, i) => (
          <MenuItem key={i} value={i}>
            {i}
          </MenuItem>
        ))}
      </Select>
    </FormControl>
  )
})
