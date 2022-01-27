import React from 'react'
import { useSelector, useDispatch } from 'react-redux'
import { FileSelect } from 'components/FlowChart/FlowChartNode/FileSelect'
import { SelectedItemIdContext } from '../VisualizeItemEditor'
import {
  selectCsvItemTranspose,
  selectVisualizeDataFilePath,
} from 'store/slice/VisualizeItem/VisualizeItemSelectors'
import {
  setCsvItemTranspose,
  setDisplayDataPath,
} from 'store/slice/VisualizeItem/VisualizeItemSlice'
import { useFileUploader } from 'store/slice/FileUploader/FileUploaderHook'
import { FILE_TYPE_SET } from 'store/slice/InputNode/InputNodeType'
import { FILE_TREE_TYPE_SET } from 'store/slice/FilesTree/FilesTreeType'
import { FormControlLabel, Switch } from '@material-ui/core'

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
    <>
      <FileSelect
        filePath={filePath ?? ''}
        onSelectFile={onSelectFile}
        onUploadFile={onUploadFileHandle}
        fileTreeType={FILE_TREE_TYPE_SET.CSV}
        selectButtonLabel="Select CSV"
      />
      <Transpose />
    </>
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
