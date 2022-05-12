import React from 'react'
import { useSelector, useDispatch } from 'react-redux'
import PlotlyChart from 'react-plotlyjs-ts'
import LinearProgress from '@mui/material/LinearProgress'
import Typography from '@mui/material/Typography'
import Box from '@mui/material/Box'
import InputLabel from '@mui/material/InputLabel'
import MenuItem from '@mui/material/MenuItem'
import FormControl from '@mui/material/FormControl'
import Select, { SelectChangeEvent } from '@mui/material/Select'

import { DisplayDataContext } from '../DataContext'
import {
  selectBarData,
  selectBarDataError,
  selectBarDataIsFulfilled,
  selectBarDataIsInitialized,
  selectBarDataIsPending,
} from 'store/slice/DisplayData/DisplayDataSelectors'
import { getBarData } from 'store/slice/DisplayData/DisplayDataActions'
import { BarData } from 'api/outputs/Outputs'
import {
  selectBarItemIndex,
  selectVisualizeItemHeight,
  selectVisualizeItemWidth,
  selectVisualizeSaveFilename,
  selectVisualizeSaveFormat,
} from 'store/slice/VisualizeItem/VisualizeItemSelectors'
import { setBarItemIndex } from 'store/slice/VisualizeItem/VisualizeItemSlice'

export const BarPlot = React.memo(() => {
  const { filePath: path } = React.useContext(DisplayDataContext)
  const dispatch = useDispatch()
  const isPending = useSelector(selectBarDataIsPending(path))
  const isInitialized = useSelector(selectBarDataIsInitialized(path))
  const error = useSelector(selectBarDataError(path))
  const isFulfilled = useSelector(selectBarDataIsFulfilled(path))
  React.useEffect(() => {
    if (!isInitialized) {
      dispatch(getBarData({ path }))
    }
  }, [dispatch, isInitialized, path])
  if (isPending) {
    return <LinearProgress />
  } else if (error != null) {
    return <Typography color="error">{error}</Typography>
  } else if (isFulfilled) {
    return <BarPlotImple />
  } else {
    return null
  }
})

const BarPlotImple = React.memo(() => {
  const { filePath: path, itemId } = React.useContext(DisplayDataContext)
  const barData = useSelector(selectBarData(path), barDataEqualityFn)
  const width = useSelector(selectVisualizeItemWidth(itemId))
  const height = useSelector(selectVisualizeItemHeight(itemId))
  const index = useSelector(selectBarItemIndex(itemId))

  const data = React.useMemo(
    () => [
      {
        x: Object.keys(barData[index]),
        y: Object.values(barData[index]),
        type: 'bar',
      },
    ],
    [barData, index],
  )

  const layout = React.useMemo(
    () => ({
      width: width,
      height: height - 120,
      margin: {
        t: 60, // top
        l: 50, // left
        b: 30, // bottom
      },
      dragmode: 'pan',
      autosize: true,
    }),
    [width, height],
  )

  const saveFileName = useSelector(selectVisualizeSaveFilename(itemId))
  const saveFormat = useSelector(selectVisualizeSaveFormat(itemId))

  const config = {
    displayModeBar: true,
    // scrollZoom: true,
    responsive: true,
    toImageButtonOptions: {
      format: saveFormat,
      filename: saveFileName,
    },
  }

  return (
    <div>
      <Box sx={{ display: 'flex' }}>
        <Box sx={{ flexGrow: 1, ml: 1 }}>
          <SelectIndex dataKeys={Object.keys(barData)} />
        </Box>
      </Box>
      <PlotlyChart data={data} layout={layout} config={config} />
    </div>
  )
})

const SelectIndex = React.memo<{
  dataKeys: string[]
}>(({ dataKeys }) => {
  const { itemId } = React.useContext(DisplayDataContext)
  const dispatch = useDispatch()
  const index = useSelector(selectBarItemIndex(itemId))

  const handleChange = (event: SelectChangeEvent<string>) => {
    dispatch(
      setBarItemIndex({
        itemId,
        index: event.target.value,
      }),
    )
  }
  return (
    <FormControl variant="standard" sx={{ m: 1, minWidth: 120 }}>
      <InputLabel>index</InputLabel>
      <Select label="smooth" value={`${index}`} onChange={handleChange}>
        {dataKeys.map((_, i) => (
          <MenuItem value={i}>{i}</MenuItem>
        ))}
      </Select>
    </FormControl>
  )
})

function barDataEqualityFn(a: BarData | undefined, b: BarData | undefined) {
  if (a != null && b != null) {
    const aArray = Object.keys(a)
    const bArray = Object.keys(b)
    return (
      a === b ||
      (aArray.length === bArray.length &&
        aArray.every((aKey, i) => {
          const bKey = bArray[i]
          return bKey === aKey // && nestEqualityFun(bValue, aValue)
        }))
    )
  } else {
    return a === undefined && b === undefined
  }
}
