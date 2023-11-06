import { memo, useContext, useEffect } from "react"
import { useSelector, useDispatch } from "react-redux"

import { LinearProgress, Typography } from "@mui/material"

import { DisplayDataContext } from "components/Workspace/Visualize/DataContext"
import { getHTMLData } from "store/slice/DisplayData/DisplayDataActions"
import {
  selectHTMLData,
  selectHTMLDataError,
  selectHTMLDataIsFulfilled,
  selectHTMLDataIsInitialized,
  selectHTMLDataIsPending,
} from "store/slice/DisplayData/DisplayDataSelectors"
import { AppDispatch } from "store/store"

export const HTMLPlot = memo(function HTMLPlot() {
  const { filePath: path } = useContext(DisplayDataContext)
  const dispatch = useDispatch<AppDispatch>()
  const isPending = useSelector(selectHTMLDataIsPending(path))
  const isInitialized = useSelector(selectHTMLDataIsInitialized(path))
  const error = useSelector(selectHTMLDataError(path))
  const isFulfilled = useSelector(selectHTMLDataIsFulfilled(path))
  useEffect(() => {
    if (!isInitialized) {
      dispatch(getHTMLData({ path }))
    }
  }, [dispatch, isInitialized, path])
  if (isPending) {
    return <LinearProgress />
  } else if (error != null) {
    return <Typography color="error">{error}</Typography>
  } else if (isFulfilled) {
    return <HTMLPlotImple />
  } else {
    return null
  }
})

const HTMLPlotImple = memo(function HTMLPlotImple() {
  const { filePath: path } = useContext(DisplayDataContext)
  const htmlData = useSelector(selectHTMLData(path))

  return (
    <div
      dangerouslySetInnerHTML={{ __html: htmlData }}
      style={{
        overflow: "scroll",
        display: "flex",
        height: "90%",
      }}
    />
  )
})
