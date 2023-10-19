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
import { getPolarData } from "store/slice/DisplayData/DisplayDataActions"
import {
  selectPolarColumns,
  selectPolarData,
  selectPolarDataError,
  selectPolarDataIsFulfilled,
  selectPolarDataIsInitialized,
  selectPolarDataIsPending,
  selectPolarIndex,
  selectPolarMeta,
} from "store/slice/DisplayData/DisplayDataSelectors"
import {
  selectPolarItemSelectedIndex,
  selectVisualizeItemHeight,
  selectVisualizeItemWidth,
} from "store/slice/VisualizeItem/VisualizeItemSelectors"
import { setPolartemItemSelectedIndex } from "store/slice/VisualizeItem/VisualizeItemSlice"
import { AppDispatch } from "store/store"

export const PolarPlot = memo(function PolarPlot() {
  const { filePath: path } = useContext(DisplayDataContext)
  const dispatch = useDispatch<AppDispatch>()
  const isPending = useSelector(selectPolarDataIsPending(path))
  const isInitialized = useSelector(selectPolarDataIsInitialized(path))
  const error = useSelector(selectPolarDataError(path))
  const isFulfilled = useSelector(selectPolarDataIsFulfilled(path))

  useEffect(() => {
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

const PolarPlotImple = memo(function PolarPlotImple() {
  const { filePath: path, itemId } = useContext(DisplayDataContext)
  const polarData = useSelector(selectPolarData(path))
  const meta = useSelector(selectPolarMeta(path))
  const columns = useSelector(selectPolarColumns(path))
  const index = useSelector(selectPolarIndex(path))
  const selectedIndex = useSelector(selectPolarItemSelectedIndex(itemId))
  const width = useSelector(selectVisualizeItemWidth(itemId))
  const height = useSelector(selectVisualizeItemHeight(itemId))

  const data = useMemo(
    () =>
      polarData != null
        ? [
            {
              type: "scatterpolar",
              mode: "lines+markae",
              theta: [...columns, columns[0]],
              r: [...polarData[selectedIndex], polarData[selectedIndex][0]],
            },
          ]
        : [],
    [polarData, columns, selectedIndex],
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
        tickvals: columns,
        title: meta?.xlabel,
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
          <PolarItemIndexSelect index={index} />
        </Box>
      </Box>
      <PlotlyChart data={data} layout={layout} />
    </>
  )
})

interface PolarItemIndexSelectProps {
  index: number[]
}

const PolarItemIndexSelect = memo(function PolarItemIndexSelect({
  index,
}: PolarItemIndexSelectProps) {
  const { itemId } = useContext(DisplayDataContext)
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
