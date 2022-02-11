import React from 'react'
import { useSelector, useDispatch } from 'react-redux'
import { styled } from '@mui/material/styles'
import Paper from '@mui/material/Paper'
import Grid from '@mui/material/Grid'

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
import { DATA_TYPE_SET } from 'store/slice/DisplayData/DisplayDataType'
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
import { FilePathSelect } from './FilePathSelect'

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

const DefaultPlot = React.memo<{
  itemId: number
}>(({ itemId }) => {
  return (
    <Grid container>
      <Grid item xs={7}>
        <StyledPaper>
          <DefaultImagePlot itemId={itemId} />
        </StyledPaper>
      </Grid>
      <Grid item xs={5}>
        <Grid>
          <StyledPaper>
            <DefaultTimeSeriesPlot itemId={itemId} />
          </StyledPaper>
          <StyledPaper>
            <DefaultHeatMapPlot itemId={itemId} />
          </StyledPaper>
        </Grid>
      </Grid>
    </Grid>
  )
})

const StyledPaper = styled(Paper)(({ theme }) => ({
  padding: theme.spacing(1),
  margin: theme.spacing(1),
  textAlign: 'center',
  color: theme.palette.text.secondary,
}))

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
