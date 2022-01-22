import React from 'react'
import { useDispatch, useSelector } from 'react-redux'

import FormControl from '@material-ui/core/FormControl'
import MenuItem from '@material-ui/core/MenuItem'
import InputLabel from '@material-ui/core/InputLabel'
import Select from '@material-ui/core/Select'
import Box from '@material-ui/core/Box'

import {
  selectSelectedVisualizeItemId,
  selectVisualizeDataFilePath,
  selectVisualizeDataNodeId,
  selectVisualizeDataType,
  selectVisualizeItemType,
} from 'store/slice/VisualizeItem/VisualizeItemSelectors'
import { VISUALIZE_ITEM_TYPE_SET } from 'store/slice/VisualizeItem/VisualizeItemType'
import {
  DATA_TYPE,
  DATA_TYPE_SET,
} from 'store/slice/DisplayData/DisplayDataType'
import { RootState } from 'store/store'
import {
  setDisplayDataPath,
  setItemType,
} from 'store/slice/VisualizeItem/VisualizeItemSlice'
import { ImageItemEditor } from './Editor/ImageItemEditor'
import { TableItemEditor } from './Editor/TableItemEditor'
import { HeatmapItemEditor } from './Editor/HeatmapItemEditor'
import { TimeSeriesItemEditor } from './Editor/TimeSeriesItemEditor'
import { RoiItemEditor } from './Editor/RoiItemEditor'
import { FilePathSelect } from './FilePathSelect'

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
  const handleChange = (event: React.ChangeEvent<{ value: unknown }>) => {
    dispatch(
      setItemType({
        itemId,
        type: event.target.value as
          | typeof VISUALIZE_ITEM_TYPE_SET.DEFAULT_SET
          | DATA_TYPE,
      }),
    )
  }
  const selectedType = useSelector((state: RootState) => {
    const itemType = selectVisualizeItemType(itemId)(state)
    if (itemType === VISUALIZE_ITEM_TYPE_SET.DEFAULT_SET) {
      return VISUALIZE_ITEM_TYPE_SET.DEFAULT_SET
    } else {
      return selectVisualizeDataType(itemId)(state)
    }
  })
  const options: (typeof VISUALIZE_ITEM_TYPE_SET.DEFAULT_SET | DATA_TYPE)[] = [
    VISUALIZE_ITEM_TYPE_SET.DEFAULT_SET,
    ...Object.values(DATA_TYPE_SET),
  ]
  return (
    <FormControl style={{ minWidth: 120, marginBottom: 12 }}>
      <InputLabel id="demo-simple-select-helper-label">item type</InputLabel>
      <Select
        value={selectedType != null ? selectedType : 'none'}
        onChange={handleChange}
      >
        {options.map((option) => (
          <MenuItem key={option} value={option}>
            {option}
          </MenuItem>
        ))}
      </Select>
    </FormControl>
  )
}

const Editor: React.FC = () => {
  const itemId = React.useContext(SelectedItemIdContext)
  const itemType = useSelector(selectVisualizeItemType(itemId))
  switch (itemType) {
    case VISUALIZE_ITEM_TYPE_SET.DEFAULT_SET:
      return <DefaultSetItemEditor />
    case VISUALIZE_ITEM_TYPE_SET.DISPLAY_DATA:
      return <DisplayDataItemEditor />
  }
}

const DefaultSetItemEditor: React.FC = () => {
  // const itemId = React.useContext(SelectedItemIdContext)
  return <div>DefaultSetItemEditor(not imple)</div>
}

const DisplayDataItemEditor: React.FC = () => {
  const itemId = React.useContext(SelectedItemIdContext)
  const dataType = useSelector(selectVisualizeDataType(itemId))
  const selectedNodeId = useSelector(selectVisualizeDataNodeId(itemId))
  const selectedFilePath = useSelector(selectVisualizeDataFilePath(itemId))
  const dispatch = useDispatch()
  const onSelect = (nodeId: string, filePath: string) =>
    dispatch(setDisplayDataPath({ itemId, nodeId, filePath }))
  return (
    <div>
      <FilePathSelect
        dataType={dataType}
        selectedNodeId={selectedNodeId}
        selectedFilePath={selectedFilePath}
        onSelect={onSelect}
      />
      <div style={{ marginTop: 8 }}>
        {dataType === DATA_TYPE_SET.IMAGE && <ImageItemEditor />}
        {/* 他のtypeのEditorも必要になったら追加する */}
        {dataType === DATA_TYPE_SET.TABLE && <TableItemEditor />}
        {dataType === DATA_TYPE_SET.HEAT_MAP && <HeatmapItemEditor />}
        {dataType === DATA_TYPE_SET.TIME_SERIES && <TimeSeriesItemEditor />}
        {dataType === DATA_TYPE_SET.ROI && <RoiItemEditor />}
      </div>
    </div>
  )
}
