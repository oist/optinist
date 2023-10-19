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

import { ScatterData } from "api/outputs/Outputs"
import { DisplayDataContext } from "components/Workspace/Visualize/DataContext"
import { getScatterData } from "store/slice/DisplayData/DisplayDataActions"
import {
  selectScatterData,
  selectScatterDataError,
  selectScatterDataIsFulfilled,
  selectScatterDataIsInitialized,
  selectScatterDataIsPending,
  selectScatterMeta,
} from "store/slice/DisplayData/DisplayDataSelectors"
import {
  selectScatterItemXIndex,
  selectScatterItemYIndex,
  selectVisualizeItemHeight,
  selectVisualizeItemWidth,
  selectVisualizeSaveFilename,
  selectVisualizeSaveFormat,
} from "store/slice/VisualizeItem/VisualizeItemSelectors"
import {
  setScatterItemXIndex,
  setScatterItemYIndex,
} from "store/slice/VisualizeItem/VisualizeItemSlice"
import { AppDispatch } from "store/store"

export const ScatterPlot = memo(function ScatterPlot() {
  const { filePath: path } = useContext(DisplayDataContext)
  const dispatch = useDispatch<AppDispatch>()
  const isPending = useSelector(selectScatterDataIsPending(path))
  const isInitialized = useSelector(selectScatterDataIsInitialized(path))
  const error = useSelector(selectScatterDataError(path))
  const isFulfilled = useSelector(selectScatterDataIsFulfilled(path))
  useEffect(() => {
    if (!isInitialized) {
      dispatch(getScatterData({ path }))
    }
  }, [dispatch, isInitialized, path])
  if (isPending) {
    return <LinearProgress />
  } else if (error != null) {
    return <Typography color="error">{error}</Typography>
  } else if (isFulfilled) {
    return <ScatterPlotImple />
  } else {
    return null
  }
})

const ScatterPlotImple = memo(function ScatterPlotImple() {
  const { filePath: path, itemId } = useContext(DisplayDataContext)

  const scatterData = useSelector(
    selectScatterData(path),
    scatterDataEqualityFn,
  )
  const meta = useSelector(selectScatterMeta(path))
  const xIndex = useSelector(selectScatterItemXIndex(itemId))
  const yIndex = useSelector(selectScatterItemYIndex(itemId))
  const width = useSelector(selectVisualizeItemWidth(itemId))
  const height = useSelector(selectVisualizeItemHeight(itemId))

  const data = useMemo(
    () => [
      {
        x: scatterData[xIndex],
        y: scatterData[yIndex],
        type: "scatter",
        mode: "markers", //'markers+text',
        text: Object.keys(scatterData[xIndex]),
        textposition: "top center",
        textfont: {
          family: "Raleway, sans-serif",
        },
        marker: {
          size: 5,
          color: Object.keys(scatterData[xIndex]),
        },
      },
    ],
    [scatterData, xIndex, yIndex],
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
        title: {
          text: meta?.xlabel ?? `x: ${xIndex}`,
          font: {
            family: "Courier New, monospace",
            size: 18,
            color: "#7f7f7f",
          },
        },
      },
      yaxis: {
        title: {
          text: meta?.ylabel ?? `y: ${yIndex}`,
          font: {
            family: "Courier New, monospace",
            size: 18,
            color: "#7f7f7f",
          },
        },
      },
    }),
    [meta, xIndex, yIndex, width, height],
  )

  const saveFileName = useSelector(selectVisualizeSaveFilename(itemId))
  const saveFormat = useSelector(selectVisualizeSaveFormat(itemId))

  const config = {
    displayModeBar: true,
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
          <XIndex dataKeys={Object.keys(scatterData)} />
          <YIndex dataKeys={Object.keys(scatterData)} />
        </Box>
      </Box>
      <PlotlyChart data={data} layout={layout} config={config} />
    </div>
  )
})

interface DataKeysProps {
  dataKeys: string[]
}

const XIndex = memo(function XIndex({ dataKeys }: DataKeysProps) {
  const { itemId } = useContext(DisplayDataContext)
  const dispatch = useDispatch()
  const xIndex = useSelector(selectScatterItemXIndex(itemId))

  const handleChange = (event: SelectChangeEvent<string>) => {
    dispatch(
      setScatterItemXIndex({
        itemId,
        xIndex: event.target.value,
      }),
    )
  }

  return (
    <FormControl variant="standard" sx={{ m: 1, minWidth: 120 }}>
      <InputLabel>xIndex</InputLabel>
      <Select label="smooth" value={xIndex} onChange={handleChange}>
        {dataKeys.map((x) => (
          <MenuItem key={x} value={x}>
            {x}
          </MenuItem>
        ))}
      </Select>
    </FormControl>
  )
})

const YIndex = memo(function YIndex({ dataKeys }: DataKeysProps) {
  const { itemId } = useContext(DisplayDataContext)
  const dispatch = useDispatch()
  const yIndex = useSelector(selectScatterItemYIndex(itemId))

  const handleChange = (event: SelectChangeEvent<string>) => {
    dispatch(
      setScatterItemYIndex({
        itemId,
        yIndex: event.target.value,
      }),
    )
  }

  return (
    <FormControl variant="standard" sx={{ m: 1, minWidth: 120 }}>
      <InputLabel>yIndex</InputLabel>
      <Select label="smooth" value={yIndex} onChange={handleChange}>
        {dataKeys.map((x) => (
          <MenuItem key={x} value={x}>
            {x}
          </MenuItem>
        ))}
      </Select>
    </FormControl>
  )
})

function scatterDataEqualityFn(
  a: ScatterData | undefined,
  b: ScatterData | undefined,
) {
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
