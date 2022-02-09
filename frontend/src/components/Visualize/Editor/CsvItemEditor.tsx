import React from 'react'
import { useSelector, useDispatch } from 'react-redux'
import { FileSelectImple } from 'components/FlowChart/FlowChartNode/FileSelect'
import { SelectedItemIdContext } from '../VisualizeItemEditor'
import {
  selectCsvItemSetColumn,
  selectCsvItemSetIndex,
  selectCsvItemTranspose,
  selectVisualizeDataFilePath,
} from 'store/slice/VisualizeItem/VisualizeItemSelectors'
import {
  setCsvItemSetColumn,
  setCsvItemSetIndex,
  setCsvItemTranspose,
  setDisplayDataPath,
} from 'store/slice/VisualizeItem/VisualizeItemSlice'
import { useFileUploader } from 'store/slice/FileUploader/FileUploaderHook'
import { FILE_TYPE_SET } from 'store/slice/InputNode/InputNodeType'
import { FILE_TREE_TYPE_SET } from 'store/slice/FilesTree/FilesTreeType'
import { FormControlLabel, Switch, TextField } from '@material-ui/core'

export const CsvItemEditor: React.FC = () => {
  const itemId = React.useContext(SelectedItemIdContext)
  const filePath = useSelector(selectVisualizeDataFilePath(itemId))
  const dispatch = useDispatch()
  const onSelectFile = (path: string) => {
    dispatch(setDisplayDataPath({ nodeId: null, filePath: path, itemId }))
  }
  const { onUploadFile } = useFileUploader(FILE_TYPE_SET.CSV)
  const onUploadFileHandle = (formData: FormData, fileName: string) => {
    onUploadFile(formData, fileName)
  }

  return (
    <div style={{ margin: '10px', padding: 10 }}>
      <FileSelectImple
        filePath={filePath ?? ''}
        onSelectFile={onSelectFile}
        onUploadFile={onUploadFileHandle}
        fileTreeType={FILE_TREE_TYPE_SET.CSV}
        selectButtonLabel="Select CSV"
      />
      <Transpose />
      <SetColumn />
      <SetIndex />
    </div>
  )
}

const Transpose: React.FC = () => {
  const itemId = React.useContext(SelectedItemIdContext)
  const transpose = useSelector(selectCsvItemTranspose(itemId))
  const dispatch = useDispatch()
  const toggleChecked = () => {
    dispatch(setCsvItemTranspose({ itemId, transpose: !transpose }))
  }
  return (
    <FormControlLabel
      control={<Switch checked={transpose} onChange={toggleChecked} />}
      label="transpose"
    />
  )
}

const SetColumn: React.FC = () => {
  const itemId = React.useContext(SelectedItemIdContext)
  const setColumn = useSelector(selectCsvItemSetColumn(itemId))

  const dispatch = useDispatch()
  const onChangeSetColumn = (event: React.ChangeEvent<HTMLInputElement>) => {
    const newValue =
      event.target.value === '' ? null : Number(event.target.value)
    if (newValue === null || newValue >= 0) {
      dispatch(setCsvItemSetColumn({ itemId, setColumn: newValue }))
    }
  }

  return (
    <FormControlLabel
      control={
        <>
          <TextField
            style={{ width: 50 }}
            type="number"
            InputLabelProps={{
              shrink: true,
            }}
            onChange={onChangeSetColumn}
            defaultValue={setColumn}
          />
          column
        </>
      }
      label=""
    />
  )
}

const SetIndex: React.FC = () => {
  const itemId = React.useContext(SelectedItemIdContext)
  const setIndex = useSelector(selectCsvItemSetIndex(itemId))
  const dispatch = useDispatch()
  const toggleChecked = () => {
    dispatch(setCsvItemSetIndex({ itemId, setIndex: !setIndex }))
  }
  return (
    <FormControlLabel
      control={<Switch checked={setIndex} onChange={toggleChecked} />}
      label="setIndex"
    />
  )
}
