import React from 'react'
import { useSelector, useDispatch } from 'react-redux'
import PlotlyChart from 'react-plotlyjs-ts'
import { DisplayDataContext } from '../DataContext'
import {
  selectPolarColumns,
  selectPolarData,
  selectPolarDataError,
  selectPolarDataIsFulfilled,
  selectPolarDataIsInitialized,
  selectPolarDataIsPending,
  selectPolarIndex,
} from 'store/slice/DisplayData/DisplayDataSelectors'
import { getPolarData } from 'store/slice/DisplayData/DisplayDataActions'
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
  selectPolarItemSelectedIndex,
  selectVisualizeItemHeight,
  selectVisualizeItemWidth,
} from 'store/slice/VisualizeItem/VisualizeItemSelectors'
import { setPolartemItemSelectedIndex } from 'store/slice/VisualizeItem/VisualizeItemSlice'

export const PolarPlot = React.memo(() => {
  const { filePath: path } = React.useContext(DisplayDataContext)
  const dispatch = useDispatch()
  const isPending = useSelector(selectPolarDataIsPending(path))
  const isInitialized = useSelector(selectPolarDataIsInitialized(path))
  const error = useSelector(selectPolarDataError(path))
  const isFulfilled = useSelector(selectPolarDataIsFulfilled(path))

  React.useEffect(() => {
    if (!isInitialized) {
      dispatch(getPolarData({ path }))
    }
  }, [dispatch, isInitialized, path])

  if (isPending) {
    return <LinearProgress />
  } else if (error != null) {
    return <Typography color="error">{error}</Typography>
  } else if (isFulfilled) {
    return <PolarPlotImple />
  } else {
    return null
  }
})

const PolarPlotImple = React.memo(() => {
  const { filePath: path, itemId } = React.useContext(DisplayDataContext)
  const polarData = useSelector(selectPolarData(path))
  const columns = useSelector(selectPolarColumns(path))
  const index = useSelector(selectPolarIndex(path))
  const selectedIndex = useSelector(selectPolarItemSelectedIndex(itemId))
  const width = useSelector(selectVisualizeItemWidth(itemId))
  const height = useSelector(selectVisualizeItemHeight(itemId))

  const data = React.useMemo(
    () =>
      polarData != null
        ? [
            {
              type: 'scatterpolar',
              mode: 'lines+markae',
              theta: [...columns, columns[0]],
              r: [...polarData[selectedIndex], polarData[selectedIndex][0]],
            },
          ]
        : [],
    [polarData, columns, selectedIndex],
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
          <PolarItemIndexSelect index={index} />
        </Box>
      </Box>
      <PlotlyChart data={data} layout={layout} />
    </>
  )
})

const PolarItemIndexSelect = React.memo<{ index: number[] }>(({ index }) => {
  const { itemId } = React.useContext(DisplayDataContext)
  const dispatch = useDispatch()
  const selectedIndex = useSelector(selectPolarItemSelectedIndex(itemId))

  const handleChange = (event: SelectChangeEvent<number>) => {
    dispatch(
      setPolartemItemSelectedIndex({
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
