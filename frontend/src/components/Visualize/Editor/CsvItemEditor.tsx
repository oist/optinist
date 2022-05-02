import React from 'react'
import { useSelector, useDispatch } from 'react-redux'
import { FileSelectImple } from 'components/FlowChart/FlowChartNode/FileSelect'
import { SelectedItemIdContext } from '../VisualizeItemEditor'
import {
  selectCsvItemSetHeader,
  selectCsvItemSetIndex,
  selectCsvItemTranspose,
  selectDisplayDataIsSingle,
  selectVisualizeDataFilePath,
} from 'store/slice/VisualizeItem/VisualizeItemSelectors'
import {
  setCsvItemSetHeader,
  setCsvItemSetIndex,
  setCsvItemTranspose,
} from 'store/slice/VisualizeItem/VisualizeItemSlice'
import { useFileUploader } from 'store/slice/FileUploader/FileUploaderHook'
import { FILE_TYPE_SET } from 'store/slice/InputNode/InputNodeType'
import { FILE_TREE_TYPE_SET } from 'store/slice/FilesTree/FilesTreeType'
import { FormControlLabel, Switch, TextField } from '@mui/material'
import { setNewDisplayDataPath } from 'store/slice/VisualizeItem/VisualizeItemActions'
import { DATA_TYPE_SET } from 'store/slice/DisplayData/DisplayDataType'

export const CsvItemEditor: React.FC = () => {
  const itemId = React.useContext(SelectedItemIdContext)
  const filePath = useSelector(selectVisualizeDataFilePath(itemId))
  const dispatch = useDispatch()
  const isSingleData = useSelector(selectDisplayDataIsSingle(itemId))
  const onSelectFile = (newPath: string) => {
    const basePayload = {
      itemId,
      nodeId: null,
      filePath: newPath,
    }
    dispatch(
      setNewDisplayDataPath(
        isSingleData && filePath != null
          ? {
              ...basePayload,
              deleteData: true,
              prevDataType: DATA_TYPE_SET.CSV,
              prevFilePath: filePath,
            }
          : {
              ...basePayload,
              deleteData: false,
            },
      ),
    )
  }
  const { onUploadFile } = useFileUploader({ fileType: FILE_TYPE_SET.CSV })
  const onUploadFileHandle = (formData: FormData, fileName: string) => {
    onUploadFile(formData, fileName)
  }

  return (
    <div style={{ margin: '10px', padding: 10 }}>
      <FileSelectImple
        filePath={filePath ?? ''}
        onSelectFile={(path) => !Array.isArray(path) && onSelectFile(path)}
        onUploadFile={onUploadFileHandle}
        fileTreeType={FILE_TREE_TYPE_SET.CSV}
        selectButtonLabel="Select CSV"
      />
      <Transpose />
      <SetHeader />
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

const SetHeader: React.FC = () => {
  const itemId = React.useContext(SelectedItemIdContext)
  const setHeader = useSelector(selectCsvItemSetHeader(itemId))

  const dispatch = useDispatch()
  const onChangeSetHeader = (event: React.ChangeEvent<HTMLInputElement>) => {
    const newValue =
      event.target.value === '' ? null : Number(event.target.value)
    if (newValue === null || newValue >= 0) {
      dispatch(setCsvItemSetHeader({ itemId, setHeader: newValue }))
    }
  }

  return (
    <>
      <TextField
        label="header"
        sx={{
          width: 100,
          margin: (theme) => theme.spacing(0, 1, 0, 1),
        }}
        type="number"
        InputLabelProps={{
          shrink: true,
        }}
        onChange={onChangeSetHeader}
        value={setHeader}
      />
    </>
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
