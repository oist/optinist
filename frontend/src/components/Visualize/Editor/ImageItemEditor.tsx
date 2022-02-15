import React, { useState } from 'react'
import { useDispatch, useSelector } from 'react-redux'
import Switch from '@mui/material/Switch'
import FormControlLabel from '@mui/material/FormControlLabel'
import Select from '@mui/material/Select'
import MenuItem from '@mui/material/MenuItem'

import {
  selectImageItemShowGrid,
  selectImageItemShowLine,
  selectImageItemShowticklabels,
  selectImageItemZsmooth,
  selectImageItemShowScale,
  selectVisualizeDataFilePath,
  selectImageItemColors,
  selectRoiItemNodeId,
  selectRoiItemFilePath,
} from 'store/slice/VisualizeItem/VisualizeItemSelectors'
import { SelectedItemIdContext } from '../VisualizeItemEditor'

import {
  setImageItemShowGrid,
  setImageItemShowLine,
  setImageItemShowticklabels,
  setImageItemZsmooth,
  setImageItemShowScale,
  setDisplayDataPath,
  setImageItemStartIndex,
  setImageItemEndIndex,
  setImageItemColors,
  setRoiItemFilePath,
  resetImageActiveIndex,
} from 'store/slice/VisualizeItem/VisualizeItemSlice'

import 'react-linear-gradient-picker/dist/index.css'
import { FileSelectImple } from 'components/FlowChart/FlowChartNode/FileSelect'
import { useFileUploader } from 'store/slice/FileUploader/FileUploaderHook'
import { FILE_TYPE_SET } from 'store/slice/InputNode/InputNodeType'
import { FILE_TREE_TYPE_SET } from 'store/slice/FilesTree/FilesTreeType'
import { TextField } from '@mui/material'
import { GradientColorPicker } from './GradientColorPicker'
import { ColorType } from 'store/slice/VisualizeItem/VisualizeItemType'
import { FilePathSelect } from '../FilePathSelect'
import { DATA_TYPE_SET } from 'store/slice/DisplayData/DisplayDataType'
import Button from '@mui/material/Button'
import { getImageData } from 'store/slice/DisplayData/DisplayDataActions'

export const ImageItemEditor: React.FC = () => {
  const itemId = React.useContext(SelectedItemIdContext)
  const dispatch = useDispatch()
  const onSelectImageFile = (path: string) => {
    dispatch(setDisplayDataPath({ nodeId: null, filePath: path, itemId }))
    dispatch(resetImageActiveIndex({ itemId }))
  }
  const filePath = useSelector(selectVisualizeDataFilePath(itemId))

  const { onUploadFile } = useFileUploader(FILE_TYPE_SET.IMAGE)
  const onUploadFileHandle = (formData: FormData, fileName: string) => {
    onUploadFile(formData, fileName)
  }

  const colors = useSelector(selectImageItemColors(itemId))
  const dispathSetColor = (colorCode: ColorType[]) => {
    dispatch(setImageItemColors({ itemId, colors: colorCode }))
  }

  const roiItemNodeId = useSelector(selectRoiItemNodeId(itemId))
  const roiItemFilePath = useSelector(selectRoiItemFilePath(itemId))
  const onSelectRoiFilePath = (nodeId: string, filePath: string) => {
    dispatch(setRoiItemFilePath({ itemId, nodeId, filePath }))
  }
  return (
    <div style={{ margin: '10px', padding: 10 }}>
      <FileSelectImple
        filePath={filePath ?? ''}
        onSelectFile={onSelectImageFile}
        onUploadFile={onUploadFileHandle}
        fileTreeType={FILE_TREE_TYPE_SET.IMAGE}
        selectButtonLabel="Select Image"
      />
      <FilePathSelect
        selectedFilePath={roiItemFilePath}
        selectedNodeId={roiItemNodeId}
        onSelect={onSelectRoiFilePath}
        dataType={DATA_TYPE_SET.ROI}
        label={'Select Roi'}
      />
      <StartEndIndex />
      <Showticklabels />
      <ShowLine />
      <ShowGrid />
      <ShowScale />
      <Zsmooth />
      <GradientColorPicker colors={colors} dispatchSetColor={dispathSetColor} />
    </div>
  )
}

const Showticklabels: React.FC = () => {
  const itemId = React.useContext(SelectedItemIdContext)
  const showticklabels = useSelector(selectImageItemShowticklabels(itemId))
  const dispatch = useDispatch()
  const toggleChecked = () => {
    dispatch(
      setImageItemShowticklabels({ itemId, showticklabels: !showticklabels }),
    )
  }
  return (
    <FormControlLabel
      control={<Switch checked={showticklabels} onChange={toggleChecked} />}
      label="Showticklabels"
    />
  )
}

const ShowLine: React.FC = () => {
  const itemId = React.useContext(SelectedItemIdContext)
  const showline = useSelector(selectImageItemShowLine(itemId))
  const dispatch = useDispatch()
  const toggleChecked = () => {
    dispatch(setImageItemShowLine({ itemId, showline: !showline }))
  }
  return (
    <FormControlLabel
      control={<Switch checked={showline} onChange={toggleChecked} />}
      label="ShowLine"
    />
  )
}

const ShowGrid: React.FC = () => {
  const itemId = React.useContext(SelectedItemIdContext)
  const showgrid = useSelector(selectImageItemShowGrid(itemId))
  const dispatch = useDispatch()
  const toggleChecked = () => {
    dispatch(setImageItemShowGrid({ itemId, showgrid: !showgrid }))
  }
  return (
    <FormControlLabel
      control={<Switch checked={showgrid} onChange={toggleChecked} />}
      label="ShowGrid"
    />
  )
}

const ShowScale: React.FC = () => {
  const itemId = React.useContext(SelectedItemIdContext)
  const showscale = useSelector(selectImageItemShowScale(itemId))
  const dispatch = useDispatch()
  const toggleChecked = () => {
    dispatch(setImageItemShowScale({ itemId, showscale: !showscale }))
  }
  return (
    <FormControlLabel
      control={<Switch checked={showscale} onChange={toggleChecked} />}
      label="ShowScale"
    />
  )
}

const Zsmooth: React.FC = () => {
  const itemId = React.useContext(SelectedItemIdContext)
  const zsmooth = useSelector(selectImageItemZsmooth(itemId))
  const dispatch = useDispatch()
  const [value, setValue] = React.useState(zsmooth)
  const handleChange = (event: any) => {
    setValue(event.target.value as string)
    dispatch(setImageItemZsmooth({ itemId, zsmooth: event.target.value }))
  }
  return (
    <FormControlLabel
      control={
        <Select
          labelId="demo-simple-select-label"
          id="demo-simple-select"
          value={value}
          onChange={handleChange}
        >
          <MenuItem value={'best'}>best</MenuItem>
          <MenuItem value={'fast'}>fast</MenuItem>
          <MenuItem value={'false'}>False</MenuItem>
        </Select>
      }
      label="smooth"
    />
  )
}

const StartEndIndex: React.FC = () => {
  const itemId = React.useContext(SelectedItemIdContext)
  const [startIndex, setStartIndex] = useState(1)
  const [endIndex, setEndIndex] = useState(10)
  const [inputError, setInputError] = useState(false)
  const dispatch = useDispatch()
  const onStartChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    const newValue = event.target.value === '' ? '' : Number(event.target.value)
    if (typeof newValue === 'number') {
      setStartIndex(newValue)
    }
  }
  const onEndChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    const newValue = event.target.value === '' ? '' : Number(event.target.value)
    if (typeof newValue === 'number') {
      setEndIndex(newValue)
    }
  }

  const filePath = useSelector(selectVisualizeDataFilePath(itemId))
  const onClickButton = () => {
    if (startIndex > 0) {
      setInputError(false)
      dispatch(resetImageActiveIndex({ itemId }))
      dispatch(setImageItemStartIndex({ itemId, startIndex: startIndex }))
      dispatch(setImageItemEndIndex({ itemId, endIndex: endIndex }))
      if (filePath !== null) {
        dispatch(
          getImageData({
            path: filePath,
            startIndex: startIndex ?? 1,
            endIndex: endIndex ?? 10,
          }),
        )
      }
    } else {
      setInputError(true)
    }
  }

  return (
    <FormControlLabel
      control={
        <>
          <TextField
            error={inputError}
            type="number"
            InputLabelProps={{
              shrink: true,
            }}
            onChange={onStartChange}
            defaultValue={startIndex}
            helperText={inputError ? 'index > 0' : ''}
          />
          ~
          <TextField
            type="number"
            InputLabelProps={{
              shrink: true,
            }}
            onChange={onEndChange}
            defaultValue={endIndex}
          />
          <Button
            size="small"
            className="ctrl_btn"
            variant="contained"
            onClick={onClickButton}
          >
            load
          </Button>
        </>
      }
      label=""
    />
  )
}
