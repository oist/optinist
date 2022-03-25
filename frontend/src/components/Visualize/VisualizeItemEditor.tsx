import React, { useState, useEffect } from 'react'
import { useDispatch, useSelector } from 'react-redux'

import FormControl from '@mui/material/FormControl'
import Box from '@mui/material/Box'
import {
  Accordion,
  AccordionDetails,
  AccordionSummary,
  FormControlLabel,
  Switch,
  Typography,
} from '@mui/material'
import ExpandMoreIcon from '@mui/icons-material/ExpandMore'

import {
  selectSelectedVisualizeItemId,
  selectImageItemFilePath,
  selectVisualizeDataNodeId,
  selectVisualizeDataType,
  selectVisualizeItemType,
  selectVisualizeItemTypeIsMultiPlot,
} from 'store/slice/VisualizeItem/VisualizeItemSelectors'
import { VISUALIZE_ITEM_TYPE_SET } from 'store/slice/VisualizeItem/VisualizeItemType'
import {
  DATA_TYPE,
  DATA_TYPE_SET,
} from 'store/slice/DisplayData/DisplayDataType'
import {
  resetImageActiveIndex,
  setDisplayDataPath,
  toggleItemTypeMultiPlot,
} from 'store/slice/VisualizeItem/VisualizeItemSlice'
import { ImageItemEditor } from './Editor/ImageItemEditor'
import { CsvItemEditor } from './Editor/CsvItemEditor'
import { HeatmapItemEditor } from './Editor/HeatmapItemEditor'
import { TimeSeriesItemEditor } from './Editor/TimeSeriesItemEditor'
import { RoiItemEditor } from './Editor/RoiItemEditor'
import { FilePathSelect } from './FilePathSelect'
import { ScatterItemEditor } from './Editor/ScatterItemEditor'
import { BarItemEditor } from './Editor/BarItemEditor'
import { deleteDisplayItem } from 'store/slice/DisplayData/DisplayDataSlice'

export const VisualizeItemEditor = () => {
  const selectedItemId = useSelector(selectSelectedVisualizeItemId)
  return (
    <>
      {selectedItemId != null ? (
        <SelectedItemIdContext.Provider value={selectedItemId}>
          <Box m={1}>
            <ItemTypeSelect />
            <Editor />
          </Box>
        </SelectedItemIdContext.Provider>
      ) : (
        'Please select item...'
      )}
      <br />
    </>
  )
}

export const SelectedItemIdContext = React.createContext<number>(NaN)

const ItemTypeSelect: React.FC = () => {
  const itemId = React.useContext(SelectedItemIdContext)
  const dispatch = useDispatch()
  const isMultiPlot = useSelector(selectVisualizeItemTypeIsMultiPlot(itemId))
  const onChageToggle = () => {
    dispatch(toggleItemTypeMultiPlot(itemId))
  }
  return (
    <FormControl style={{ minWidth: 120, marginBottom: 12 }}>
      <FormControlLabel
        control={<Switch checked={isMultiPlot} onChange={onChageToggle} />}
        label="Multi plot"
      />
    </FormControl>
  )
}

const Editor: React.FC = () => {
  const itemId = React.useContext(SelectedItemIdContext)
  const itemType = useSelector(selectVisualizeItemType(itemId))
  switch (itemType) {
    case VISUALIZE_ITEM_TYPE_SET.MULTI_PLOT:
      return <MultiPlotItemEditor />
    case VISUALIZE_ITEM_TYPE_SET.DISPLAY_DATA:
      return <DisplayDataItemEditor />
  }
}

const MultiPlotItemEditor: React.FC = () => {
  return (
    <div>
      <Accordion>
        <AccordionSummary
          expandIcon={<ExpandMoreIcon />}
          aria-controls="panel1a-content"
          id="panel1a-header"
        >
          <Typography>ImageEditor</Typography>
        </AccordionSummary>
        <AccordionDetails>
          <ImageItemEditor />
        </AccordionDetails>
      </Accordion>

      <Accordion>
        <AccordionSummary
          expandIcon={<ExpandMoreIcon />}
          aria-controls="panel2a-content"
          id="panel2a-header"
        >
          <Typography>TimeSeriesEditor</Typography>
        </AccordionSummary>
        <AccordionDetails>
          <TimeSeriesItemEditor />
        </AccordionDetails>
      </Accordion>

      <Accordion>
        <AccordionSummary
          expandIcon={<ExpandMoreIcon />}
          aria-controls="panel3a-content"
          id="panel3a-header"
        >
          <Typography>HeatmapEditor</Typography>
        </AccordionSummary>
        <AccordionDetails>
          <HeatmapItemEditor />
        </AccordionDetails>
      </Accordion>
    </div>
  )
}

const DisplayDataItemEditor: React.FC = () => {
  const itemId = React.useContext(SelectedItemIdContext)
  const dataType = useSelector(selectVisualizeDataType(itemId))
  const selectedNodeId = useSelector(selectVisualizeDataNodeId(itemId))
  const selectedFilePath = useSelector(selectImageItemFilePath(itemId))

  const [prevItem, setPrevItem] = useState<{
    dataType: DATA_TYPE
    filePath: string | null
  }>({
    dataType: 'image',
    filePath: null,
  })

  useEffect(() => {
    setPrevItem({ dataType, filePath: selectedFilePath })
  }, [selectedFilePath, dataType])

  const dispatch = useDispatch()
  const onSelect = (nodeId: string, filePath: string, dataType: DATA_TYPE) => {
    dispatch(setDisplayDataPath({ itemId, nodeId, filePath, dataType }))
    dispatch(resetImageActiveIndex({ itemId }))
    dispatch(deleteDisplayItem(prevItem))
  }

  return (
    <>
      <FilePathSelect
        selectedNodeId={selectedNodeId}
        selectedFilePath={selectedFilePath}
        onSelect={onSelect}
      />
      <div style={{ marginTop: 8 }}>
        <DisplayEditor dataType={dataType} />
      </div>
    </>
  )
}

const DisplayEditor: React.FC<{
  dataType: DATA_TYPE | null
}> = ({ dataType }) => {
  /* 他のtypeのEditorも必要になったら追加する */
  switch (dataType) {
    case DATA_TYPE_SET.IMAGE:
      return <ImageItemEditor />
    case DATA_TYPE_SET.CSV:
      return <CsvItemEditor />
    case DATA_TYPE_SET.HEAT_MAP:
      return <HeatmapItemEditor />
    case DATA_TYPE_SET.TIME_SERIES:
      return <TimeSeriesItemEditor />
    case DATA_TYPE_SET.ROI:
      return <RoiItemEditor />
    case DATA_TYPE_SET.SCATTER:
      return <ScatterItemEditor />
    case DATA_TYPE_SET.BAR:
      return <BarItemEditor />
    case DATA_TYPE_SET.HTML:
      return <div>html editor</div>
    default:
      return null
  }
}
