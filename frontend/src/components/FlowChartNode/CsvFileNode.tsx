import React, { CSSProperties } from 'react'
import { useDispatch, useSelector } from 'react-redux'
import { Handle, Position, NodeProps } from 'react-flow-renderer'
import { FlexLayoutModelContext } from 'App'
import { useTabAction } from 'FlexLayoutHook'
import { OUTPUT_TABSET_ID } from 'const/flexlayout'
import { FileSelect } from './FileSelect'
import { alpha, useTheme } from '@material-ui/core'
import { selectCsvFile } from 'store/slice/FileData/FileData'
import { uploadCsvFile } from 'store/slice/FileData/FileDataAction'
import {
  csvIsUploadingByIdSelector,
  csvUploadingProgressSelector,
} from 'store/slice/FileData/FileDataSelector'
import { FILE_TYPE_SET } from 'store/slice/FilesTree/FilesTreeType'
import { LinearProgressWithLabel } from './LinerProgressWithLabel'

const sourceHandleStyle: CSSProperties = {
  width: 8,
  height: 15,
  top: 15,
  border: '1px solid',
  borderColor: '#555',
  borderRadius: 0,
}

export const CsvFileNode = React.memo<NodeProps>(({ id: nodeId, selected }) => {
  const model = React.useContext(FlexLayoutModelContext)
  const actionForImageTab = useTabAction(nodeId)
  const onClick = () => {
    if (actionForImageTab != null) {
      model.doAction(actionForImageTab('output', OUTPUT_TABSET_ID))
    }
  }
  const theme = useTheme()
  const csvIsUploading = useSelector(csvIsUploadingByIdSelector(nodeId))
  const uploadProgress = useSelector(csvUploadingProgressSelector(nodeId))
  return (
    <div
      style={{
        height: '100%',
        background: selected
          ? alpha(theme.palette.primary.light, 0.1)
          : undefined,
      }}
      onClick={onClick}
    >
      {csvIsUploading && uploadProgress != null && (
        <div style={{ marginLeft: 2, marginRight: 2 }}>
          <LinearProgressWithLabel value={uploadProgress} />
        </div>
      )}
      <CsvFileSelect nodeId={nodeId} />
      <Handle
        type="source"
        position={Position.Right}
        id={`csv-${nodeId}`}
        style={sourceHandleStyle}
      />
    </div>
  )
})

const CsvFileSelect = React.memo<{ nodeId: string }>(({ nodeId }) => {
  const dispatch = useDispatch()
  const model = React.useContext(FlexLayoutModelContext)
  const actionForImageTab = useTabAction(nodeId)
  const onUploadFile = (formData: FormData, fileName: string) => {
    dispatch(
      uploadCsvFile({
        nodeId,
        fileName,
        formData,
      }),
    )
    if (actionForImageTab != null) {
      model.doAction(actionForImageTab('output', OUTPUT_TABSET_ID))
    }
  }
  const onSelectFile = (path: string) => {
    dispatch(selectCsvFile({ nodeId, path }))
    if (actionForImageTab != null) {
      model.doAction(actionForImageTab('output', OUTPUT_TABSET_ID))
    }
  }
  return (
    <FileSelect
      nodeId={nodeId}
      onSelectFile={onSelectFile}
      onUploadFile={onUploadFile}
      fileType={FILE_TYPE_SET.CSV}
      selectButtonLabel="csvを選択"
    />
  )
})
