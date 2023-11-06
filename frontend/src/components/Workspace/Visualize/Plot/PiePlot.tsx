import { memo, useContext, useEffect, useMemo } from "react"
import PlotlyChart from "react-plotlyjs-ts"
import { useSelector, useDispatch } from "react-redux"

import { LinearProgress, Typography } from "@mui/material"

import { DisplayDataContext } from "components/Workspace/Visualize/DataContext"
import { getPieData } from "store/slice/DisplayData/DisplayDataActions"
import {
  selectPieColumns,
  selectPieData,
  selectPieDataError,
  selectPieDataIsFulfilled,
  selectPieDataIsInitialized,
  selectPieDataIsPending,
  selectPieMeta,
} from "store/slice/DisplayData/DisplayDataSelectors"
import {
  selectVisualizeItemHeight,
  selectVisualizeItemWidth,
} from "store/slice/VisualizeItem/VisualizeItemSelectors"
import { AppDispatch } from "store/store"

export const PiePlot = memo(function PiePlot() {
  const { filePath: path } = useContext(DisplayDataContext)
  const dispatch = useDispatch<AppDispatch>()
  const isPending = useSelector(selectPieDataIsPending(path))
  const isInitialized = useSelector(selectPieDataIsInitialized(path))
  const error = useSelector(selectPieDataError(path))
  const isFulfilled = useSelector(selectPieDataIsFulfilled(path))

  useEffect(() => {
    if (!isInitialized) {
      dispatch(getPieData({ path }))
    }
  }, [dispatch, isInitialized, path])

  if (isPending) {
    return <LinearProgress />
  } else if (error != null) {
    return <Typography color="error">{error}</Typography>
  } else if (isFulfilled) {
    return <PiePlotImple />
  } else {
    return null
  }
})

const PiePlotImple = memo(function PiePlotImple() {
  const { filePath: path, itemId } = useContext(DisplayDataContext)
  const pieData = useSelector(selectPieData(path))
  const meta = useSelector(selectPieMeta(path))
  const columns = useSelector(selectPieColumns(path))
  const width = useSelector(selectVisualizeItemWidth(itemId))
  const height = useSelector(selectVisualizeItemHeight(itemId))

  const data = useMemo(
    () =>
      pieData != null
        ? [
            {
              values: pieData[0],
              labels: columns,
              type: "pie",
              sort: false,
              direction: "clockwise",
            },
          ]
        : [],
    [pieData, columns],
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
        t: 60, // top
        l: 50, // left
        b: 30, // bottom
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

  return <PlotlyChart data={data} layout={layout} />
})
