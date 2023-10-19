import { memo, useContext, useEffect, useMemo } from "react"
import PlotlyChart from "react-plotlyjs-ts"
import { useSelector, useDispatch } from "react-redux"

import Box from "@mui/material/Box"
import FormControl from "@mui/material/FormControl"
import InputLabel from "@mui/material/InputLabel"
import LinearProgress from "@mui/material/LinearProgress"
import MenuItem from "@mui/material/MenuItem"
import Select, { SelectChangeEvent } from "@mui/material/Select"
import Typography from "@mui/material/Typography"

import { BarData } from "api/outputs/Outputs"
import { DisplayDataContext } from "components/Workspace/Visualize/DataContext"
import { getBarData } from "store/slice/DisplayData/DisplayDataActions"
import {
  selectBarData,
  selectBarDataError,
  selectBarDataIsFulfilled,
  selectBarDataIsInitialized,
  selectBarDataIsPending,
  selectBarIndex,
  selectBarMeta,
} from "store/slice/DisplayData/DisplayDataSelectors"
import {
  selectBarItemIndex,
  selectVisualizeItemHeight,
  selectVisualizeItemWidth,
  selectVisualizeSaveFilename,
  selectVisualizeSaveFormat,
} from "store/slice/VisualizeItem/VisualizeItemSelectors"
import { setBarItemIndex } from "store/slice/VisualizeItem/VisualizeItemSlice"
import { AppDispatch } from "store/store"

export const BarPlot = memo(function BarPlot() {
  const { filePath: path } = useContext(DisplayDataContext)
  const dispatch = useDispatch<AppDispatch>()
  const isPending = useSelector(selectBarDataIsPending(path))
  const isInitialized = useSelector(selectBarDataIsInitialized(path))
  const error = useSelector(selectBarDataError(path))
  const isFulfilled = useSelector(selectBarDataIsFulfilled(path))
  useEffect(() => {
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

const BarPlotImple = memo(function BarPlotImple() {
  const { filePath: path, itemId } = useContext(DisplayDataContext)
  const barData = useSelector(selectBarData(path), barDataEqualityFn)
  const meta = useSelector(selectBarMeta(path))
  const width = useSelector(selectVisualizeItemWidth(itemId))
  const height = useSelector(selectVisualizeItemHeight(itemId))
  const index = useSelector(selectBarItemIndex(itemId))
  const dataKeys = useSelector(selectBarIndex(path))

  const data = useMemo(
    () => [
      {
        x: Object.keys(barData[index]),
        y: Object.values(barData[index]),
        type: "bar",
      },
    ],
    [barData, index],
  )

  const layout = useMemo(
    () => ({
      title: {
        text: meta?.title,
        x: 0.1,
      },
      width: width,
      height: height - 120,
      margin: {
        t: 50, // top
        l: 50, // left
        b: 40, // bottom
      },
      dragmode: "pan",
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
      <Box sx={{ display: "flex" }}>
        <Box sx={{ flexGrow: 1, ml: 1 }}>
          <SelectIndex dataKeys={dataKeys} />
        </Box>
      </Box>
      <PlotlyChart data={data} layout={layout} config={config} />
    </div>
  )
})

interface SelectIndexProps {
  dataKeys: string[]
}

const SelectIndex = memo(function SelectIndex({ dataKeys }: SelectIndexProps) {
  const { itemId } = useContext(DisplayDataContext)
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
        {dataKeys.map((v, i) => (
          <MenuItem key={i} value={i}>
            {v}
          </MenuItem>
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
