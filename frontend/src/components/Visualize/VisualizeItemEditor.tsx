import React from 'react'
import { useSelector } from 'react-redux'

import Box from '@mui/material/Box'

import {
  selectSelectedVisualizeItemId,
  selectVisualizeDataType,
} from 'store/slice/VisualizeItem/VisualizeItemSelectors'
import {
  DATA_TYPE,
  DATA_TYPE_SET,
} from 'store/slice/DisplayData/DisplayDataType'
import { ImageItemEditor } from './Editor/ImageItemEditor'
import { CsvItemEditor } from './Editor/CsvItemEditor'
import { HeatmapItemEditor } from './Editor/HeatmapItemEditor'
import { TimeSeriesItemEditor } from './Editor/TimeSeriesItemEditor'
import { RoiItemEditor } from './Editor/RoiItemEditor'
import { ScatterItemEditor } from './Editor/ScatterItemEditor'
import { BarItemEditor } from './Editor/BarItemEditor'

export const VisualizeItemEditor = () => {
  const selectedItemId = useSelector(selectSelectedVisualizeItemId)
  return (
    <>
      {selectedItemId != null ? (
        <SelectedItemIdContext.Provider value={selectedItemId}>
          <Box m={1}>
            <DisplayDataItemEditor />
          </Box>
        </SelectedItemIdContext.Provider>
      ) : (
        'Please select item...'
      )}
    </>
  )
}

export const SelectedItemIdContext = React.createContext<number>(NaN)

const DisplayDataItemEditor: React.FC = () => {
  const itemId = React.useContext(SelectedItemIdContext)
  const dataType = useSelector(selectVisualizeDataType(itemId))
  return (
    <div style={{ marginTop: 8 }}>
      <DisplayEditor dataType={dataType} />
    </div>
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
