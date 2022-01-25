import React from 'react'
import { useSelector, useDispatch } from 'react-redux'
import { FilePathSelect } from './FilePathSelect'
import {
  selectDefaultSetImageItemNodeId,
  selectDefaultSetImageItemFilePath,
  selectDefaultSetTimeSeriesItemNodeId,
  selectDefaultSetTimeSeriesItemFilePath,
  selectDefaultSetHeatMapItemNodeId,
  selectDefaultSetHeatMapItemFilePath,
  selectDefaultSetRoiItemNodeId,
  selectDefaultSetRoiItemFilePath,
} from 'store/slice/VisualizeItem/VisualizeItemSelectors'
import {
  setHeatMapItemFilePath,
  setImageItemFilePath,
  setTimeSeriesItemFilePath,
  setRoiItemFilePath,
} from 'store/slice/VisualizeItem/VisualizeItemSlice'
import { ImagePlot } from './Plot/ImagePlot'
import { DisplayDataContext } from './DataContext'
import { TimeSeriesPlot } from './Plot/TimeSeriesPlot'
import { HeatMapPlot } from './Plot/HeatMapPlot'

import { makeStyles, createStyles, Theme } from '@material-ui/core/styles'
import Paper from '@material-ui/core/Paper'
import Grid from '@material-ui/core/Grid'
import { DATA_TYPE_SET } from 'store/slice/DisplayData/DisplayDataType'

export const DefaultSetItem = React.memo<{
  itemId: number
}>(({ itemId }) => {
  const dispatch = useDispatch()
  return (
    <>
      <FilePathSelect
        dataType={DATA_TYPE_SET.IMAGE}
        selectedNodeId={useSelector(selectDefaultSetImageItemNodeId(itemId))}
        selectedFilePath={useSelector(
          selectDefaultSetImageItemFilePath(itemId),
        )}
        onSelect={(nodeId, filePath) =>
          dispatch(setImageItemFilePath({ itemId, nodeId, filePath }))
        }
        label="Select image"
      />
      <FilePathSelect
        dataType={DATA_TYPE_SET.TIME_SERIES}
        selectedNodeId={useSelector(
          selectDefaultSetTimeSeriesItemNodeId(itemId),
        )}
        selectedFilePath={useSelector(
          selectDefaultSetTimeSeriesItemFilePath(itemId),
        )}
        onSelect={(nodeId, filePath) =>
          dispatch(setTimeSeriesItemFilePath({ itemId, nodeId, filePath }))
        }
        label="Select timeseries"
      />
      <FilePathSelect
        dataType={DATA_TYPE_SET.HEAT_MAP}
        selectedNodeId={useSelector(selectDefaultSetHeatMapItemNodeId(itemId))}
        selectedFilePath={useSelector(
          selectDefaultSetHeatMapItemFilePath(itemId),
        )}
        onSelect={(nodeId, filePath) =>
          dispatch(setHeatMapItemFilePath({ itemId, nodeId, filePath }))
        }
        label="Select heatmap"
      />
      <FilePathSelect
        dataType={DATA_TYPE_SET.ROI}
        selectedNodeId={useSelector(selectDefaultSetRoiItemNodeId(itemId))}
        selectedFilePath={useSelector(selectDefaultSetRoiItemFilePath(itemId))}
        onSelect={(nodeId, filePath) =>
          dispatch(setRoiItemFilePath({ itemId, nodeId, filePath }))
        }
        label="Select roi"
      />
      <DefaultPlot itemId={itemId} />
    </>
  )
})

const useStyles = makeStyles((theme: Theme) =>
  createStyles({
    root: {
      flexGrow: 1,
    },
    paper: {
      padding: theme.spacing(2),
      textAlign: 'center',
      color: theme.palette.text.secondary,
    },
  }),
)

const DefaultPlot = React.memo<{
  itemId: number
}>(({ itemId }) => {
  const classes = useStyles()
  return (
    <Grid container>
      <Grid item xs={7}>
        <Paper className={classes.paper}>
          <DefaultImagePlot itemId={itemId} />
        </Paper>
      </Grid>
      <Grid item xs={5}>
        <Grid>
          <Paper className={classes.paper}>
            <DefaultTimeSeriesPlot itemId={itemId} />
          </Paper>
          <Paper className={classes.paper}>
            <DefaultHeatMapPlot itemId={itemId} />
          </Paper>
        </Grid>
      </Grid>
    </Grid>
  )
})

const DefaultImagePlot = React.memo<{
  itemId: number
}>(({ itemId }) => {
  const filePath = useSelector(selectDefaultSetImageItemFilePath(itemId))
  const nodeId = useSelector(selectDefaultSetImageItemNodeId(itemId))
  if (filePath != null) {
    return (
      <DisplayDataContext.Provider
        value={{ nodeId, filePath, dataType: DATA_TYPE_SET.IMAGE, itemId }}
      >
        <ImagePlot />
      </DisplayDataContext.Provider>
    )
  } else {
    return <div>Please select item correctly.</div>
  }
})

const DefaultTimeSeriesPlot = React.memo<{
  itemId: number
}>(({ itemId }) => {
  const filePath = useSelector(selectDefaultSetTimeSeriesItemFilePath(itemId))
  const nodeId = useSelector(selectDefaultSetTimeSeriesItemNodeId(itemId))
  if (filePath != null) {
    return (
      <DisplayDataContext.Provider
        value={{
          nodeId,
          filePath,
          dataType: DATA_TYPE_SET.TIME_SERIES,
          itemId,
        }}
      >
        <TimeSeriesPlot />
      </DisplayDataContext.Provider>
    )
  } else {
    return <div>Please select item correctly.</div>
  }
})

const DefaultHeatMapPlot = React.memo<{
  itemId: number
}>(({ itemId }) => {
  const filePath = useSelector(selectDefaultSetHeatMapItemFilePath(itemId))
  const nodeId = useSelector(selectDefaultSetHeatMapItemNodeId(itemId))
  if (filePath != null) {
    return (
      <DisplayDataContext.Provider
        value={{ nodeId, filePath, dataType: DATA_TYPE_SET.HEAT_MAP, itemId }}
      >
        <HeatMapPlot />
      </DisplayDataContext.Provider>
    )
  } else {
    return <div>Please select item correctly.</div>
  }
})
