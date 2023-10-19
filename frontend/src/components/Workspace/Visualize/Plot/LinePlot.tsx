import { memo, useContext, useEffect, useMemo } from "react"
import PlotlyChart from "react-plotlyjs-ts"
import { useSelector, useDispatch } from "react-redux"

import {
  Box,
  FormControl,
  InputLabel,
  LinearProgress,
  MenuItem,
  Select,
  SelectChangeEvent,
  Typography,
} from "@mui/material"

import { DisplayDataContext } from "components/Workspace/Visualize/DataContext"
import { getLineData } from "store/slice/DisplayData/DisplayDataActions"
import {
  selectLineColumns,
  selectLineData,
  selectLineDataError,
  selectLineDataIsFulfilled,
  selectLineDataIsInitialized,
  selectLineDataIsPending,
  selectLineIndex,
  selectLineMeta,
} from "store/slice/DisplayData/DisplayDataSelectors"
import {
  selectLineItemSelectedIndex,
  selectVisualizeItemHeight,
  selectVisualizeItemWidth,
} from "store/slice/VisualizeItem/VisualizeItemSelectors"
import { setLineItemSelectedIndex } from "store/slice/VisualizeItem/VisualizeItemSlice"
import { AppDispatch } from "store/store"

export const LinePlot = memo(function LinePlot() {
  const { filePath: path } = useContext(DisplayDataContext)
  const dispatch = useDispatch<AppDispatch>()
  const isPending = useSelector(selectLineDataIsPending(path))
  const isInitialized = useSelector(selectLineDataIsInitialized(path))
  const error = useSelector(selectLineDataError(path))
  const isFulfilled = useSelector(selectLineDataIsFulfilled(path))

  useEffect(() => {
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

const LinePlotImple = memo(function LinePlotImple() {
  const { filePath: path, itemId } = useContext(DisplayDataContext)
  const lineData = useSelector(selectLineData(path))
  const meta = useSelector(selectLineMeta(path))
  const columns = useSelector(selectLineColumns(path))
  const index = useSelector(selectLineIndex(path))
  const selectedIndex = useSelector(selectLineItemSelectedIndex(itemId))
  const width = useSelector(selectVisualizeItemWidth(itemId))
  const height = useSelector(selectVisualizeItemHeight(itemId))

  const data = useMemo(
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

  const layout = useMemo(
    () => ({
      title: {
        text: meta?.title,
        x: 0.1,
      },
      width: width,
      height: height - 120,
      dragmode: "pan",
      margin: {
        t: 50, // top
        l: 50, // left
        b: 40, // bottom
      },
      autosize: true,
      xaxis: {
        title: meta?.xlabel,
        tickvals: columns,
      },
      yaxis: {
        title: meta?.ylabel,
      },
    }),
    [meta, width, height, columns],
  )

  return (
    <>
      <Box sx={{ display: "flex" }}>
        <Box sx={{ flexGrow: 1, ml: 1 }}>
          <LineItemIndexSelect index={index} />
        </Box>
      </Box>
      <PlotlyChart data={data} layout={layout} />
    </>
  )
})

interface LineItemIndexSelectProps {
  index: number[]
}

const LineItemIndexSelect = memo(function LineItemIndexSelect({
  index,
}: LineItemIndexSelectProps) {
  const { itemId } = useContext(DisplayDataContext)
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
