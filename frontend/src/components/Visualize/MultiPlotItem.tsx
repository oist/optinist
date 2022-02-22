import React from 'react'
import { useSelector, useDispatch } from 'react-redux'
import { styled } from '@mui/material/styles'
import Paper from '@mui/material/Paper'
import Grid from '@mui/material/Grid'

import {
  selectMultiPlotImageItemNodeId,
  selectMultiPlotImageItemFilePath,
  selectMultiPlotTimeSeriesItemNodeId,
  selectMultiPlotTimeSeriesItemFilePath,
  selectMultiPlotHeatMapItemNodeId,
  selectMultiPlotHeatMapItemFilePath,
  selectMultiPlotRoiItemNodeId,
  selectMultiPlotRoiItemFilePath,
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

export const MultiPlotItem = React.memo<{
  itemId: number
}>(({ itemId }) => {
  const dispatch = useDispatch()
  return (
    <>
      <FilePathSelect
        dataType={DATA_TYPE_SET.IMAGE}
        selectedNodeId={useSelector(selectMultiPlotImageItemNodeId(itemId))}
        selectedFilePath={useSelector(selectMultiPlotImageItemFilePath(itemId))}
        onSelect={(nodeId, filePath) =>
          dispatch(setImageItemFilePath({ itemId, nodeId, filePath }))
        }
        label="Select image"
      />
      <FilePathSelect
        dataType={DATA_TYPE_SET.TIME_SERIES}
        selectedNodeId={useSelector(
          selectMultiPlotTimeSeriesItemNodeId(itemId),
        )}
        selectedFilePath={useSelector(
          selectMultiPlotTimeSeriesItemFilePath(itemId),
        )}
        onSelect={(nodeId, filePath) =>
          dispatch(setTimeSeriesItemFilePath({ itemId, nodeId, filePath }))
        }
        label="Select timeseries"
      />
      <FilePathSelect
        dataType={DATA_TYPE_SET.HEAT_MAP}
        selectedNodeId={useSelector(selectMultiPlotHeatMapItemNodeId(itemId))}
        selectedFilePath={useSelector(
          selectMultiPlotHeatMapItemFilePath(itemId),
        )}
        onSelect={(nodeId, filePath) =>
          dispatch(setHeatMapItemFilePath({ itemId, nodeId, filePath }))
        }
        label="Select heatmap"
      />
      <FilePathSelect
        dataType={DATA_TYPE_SET.ROI}
        selectedNodeId={useSelector(selectMultiPlotRoiItemNodeId(itemId))}
        selectedFilePath={useSelector(selectMultiPlotRoiItemFilePath(itemId))}
        onSelect={(nodeId, filePath) =>
          dispatch(setRoiItemFilePath({ itemId, nodeId, filePath }))
        }
        label="Select roi"
      />
      <MultiPlot itemId={itemId} />
    </>
  )
})

const MultiPlot = React.memo<{
  itemId: number
}>(({ itemId }) => {
  return (
    <Grid container>
      <Grid item xs={7}>
        <StyledPaper>
          <MultiImagePlot itemId={itemId} />
        </StyledPaper>
      </Grid>
      <Grid item xs={5}>
        <Grid>
          <StyledPaper>
            <MultiTimeSeriesPlot itemId={itemId} />
          </StyledPaper>
          <StyledPaper>
            <MultiHeatMapPlot itemId={itemId} />
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

const MultiImagePlot = React.memo<{
  itemId: number
}>(({ itemId }) => {
  const filePath = useSelector(selectMultiPlotImageItemFilePath(itemId))
  const nodeId = useSelector(selectMultiPlotImageItemNodeId(itemId))
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

const MultiTimeSeriesPlot = React.memo<{
  itemId: number
}>(({ itemId }) => {
  const filePath = useSelector(selectMultiPlotTimeSeriesItemFilePath(itemId))
  const nodeId = useSelector(selectMultiPlotTimeSeriesItemNodeId(itemId))
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

const MultiHeatMapPlot = React.memo<{
  itemId: number
}>(({ itemId }) => {
  const filePath = useSelector(selectMultiPlotHeatMapItemFilePath(itemId))
  const nodeId = useSelector(selectMultiPlotHeatMapItemNodeId(itemId))
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
