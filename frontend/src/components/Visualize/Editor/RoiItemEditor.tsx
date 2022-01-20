import React from 'react'
import { useSelector, useDispatch } from 'react-redux'
import { FileSelect } from 'components/FlowChart/FlowChartNode/FileSelect'
import { SelectedItemIdContext } from '../VisualizeItemEditor'
import {
  selectRoiItemColors,
  selectVisualizeDataFilePath,
} from 'store/slice/VisualizeItem/VisualizeItemSelectors'
import {
  setDisplayDataPath,
  setRoiItemColors,
} from 'store/slice/VisualizeItem/VisualizeItemSlice'
import { useFileUploader } from 'store/slice/FileUploader/FileUploaderHook'
import { FILE_TYPE_SET } from 'store/slice/InputNode/InputNodeType'
import { FILE_TREE_TYPE_SET } from 'store/slice/FilesTree/FilesTreeType'
import { ColorType } from 'store/slice/VisualizeItem/VisualizeItemType'
import { GradientColorPicker } from './GradientColorPicker'

export const RoiItemEditor: React.FC = () => {
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

  const colors = useSelector(selectRoiItemColors(itemId))
  const dispathSetColor = (colorCode: ColorType[]) => {
    dispatch(setRoiItemColors({ itemId, colors: colorCode }))
  }

  return (
    <>
      {/* <FileSelect
        filePath={filePath ?? ''}
        onSelectFile={onSelectFile}
        onUploadFile={onUploadFileHandle}
        fileTreeType={FILE_TREE_TYPE_SET.CSV}
        selectButtonLabel="Select CSV"
      /> */}
      {/* <GradientColorPicker /> */}
      <GradientColorPicker colors={colors} dispatchSetColor={dispathSetColor} />
    </>
  )
}
