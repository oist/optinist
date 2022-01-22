import React from 'react'
import { useSelector } from 'react-redux'
import { FilePathSelect } from './FilePathSelect'
import {
  selectDefaultSetFilePath,
  selectDefaultSetNodeId,
} from 'store/slice/VisualizeItem/VisualizeItemSelectors'
import { ImagePlot } from './Plot/ImagePlot'
import { DisplayDataContext } from './DataContext'
import { TimeSeries } from './Plot/TimeSeries'
import { HeatMap } from './Plot/HeatMap'

import { makeStyles, createStyles, Theme } from '@material-ui/core/styles'
import Paper from '@material-ui/core/Paper'
import Grid from '@material-ui/core/Grid'

export const DefaultSetItem = React.memo<{
  itemId: number
}>(({ itemId }) => {
  return (
    <>
      <FilePathSelect
        itemId={itemId}
        dataType={'image'}
        selectedNodeId={useSelector(selectDefaultSetNodeId(itemId, 'image'))}
        selectedFilePath={useSelector(
          selectDefaultSetFilePath(itemId, 'image'),
        )}
      />
      <FilePathSelect
        itemId={itemId}
        dataType={'timeSeries'}
        selectedNodeId={useSelector(
          selectDefaultSetNodeId(itemId, 'timeSeries'),
        )}
        selectedFilePath={useSelector(
          selectDefaultSetFilePath(itemId, 'timeSeries'),
        )}
      />
      <FilePathSelect
        itemId={itemId}
        dataType={'heatMap'}
        selectedNodeId={useSelector(selectDefaultSetNodeId(itemId, 'heatMap'))}
        selectedFilePath={useSelector(
          selectDefaultSetFilePath(itemId, 'heatMap'),
        )}
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
  const dataType = 'image'
  const filePath = useSelector(selectDefaultSetFilePath(itemId, dataType))
  const nodeId = useSelector(selectDefaultSetNodeId(itemId, dataType))
  if (filePath != null && dataType != null) {
    return (
      <DisplayDataContext.Provider
        value={{ nodeId, filePath, dataType, itemId }}
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
  const dataType = 'timeSeries'
  const filePath = useSelector(selectDefaultSetFilePath(itemId, dataType))
  const nodeId = useSelector(selectDefaultSetNodeId(itemId, dataType))
  if (filePath != null && dataType != null) {
    return (
      <DisplayDataContext.Provider
        value={{ nodeId, filePath, dataType, itemId }}
      >
        <TimeSeries />
      </DisplayDataContext.Provider>
    )
  } else {
    return <div>Please select item correctly.</div>
  }
})

const DefaultHeatMapPlot = React.memo<{
  itemId: number
}>(({ itemId }) => {
  const dataType = 'heatMap'
  const filePath = useSelector(selectDefaultSetFilePath(itemId, dataType))
  const nodeId = useSelector(selectDefaultSetNodeId(itemId, dataType))
  if (filePath != null && dataType != null) {
    return (
      <DisplayDataContext.Provider
        value={{ nodeId, filePath, dataType, itemId }}
      >
        <HeatMap />
      </DisplayDataContext.Provider>
    )
  } else {
    return <div>Please select item correctly.</div>
  }
})
