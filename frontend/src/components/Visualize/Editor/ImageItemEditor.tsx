import React from 'react'
import { useDispatch, useSelector } from 'react-redux'
import Switch from '@mui/material/Switch'
import FormControlLabel from '@mui/material/FormControlLabel'
import Select, { SelectChangeEvent } from '@mui/material/Select'
import InputLabel from '@mui/material/InputLabel'
import MenuItem from '@mui/material/MenuItem'
import FormControl from '@mui/material/FormControl'

import {
  selectImageItemShowGrid,
  selectImageItemShowLine,
  selectImageItemShowticklabels,
  selectImageItemZsmooth,
  selectImageItemShowScale,
  selectImageItemColors,
  selectRoiItemNodeId,
  selectRoiItemFilePath,
  selectImageItemStartIndex,
  selectImageItemEndIndex,
  selectImageItemRoiAlpha,
  selectImageItemFilePath,
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
  setImageItemRoiAlpha,
} from 'store/slice/VisualizeItem/VisualizeItemSlice'

import 'react-linear-gradient-picker/dist/index.css'
import { FileSelectImple } from 'components/FlowChart/FlowChartNode/FileSelect'
import { useFileUploader } from 'store/slice/FileUploader/FileUploaderHook'
import { FILE_TYPE_SET } from 'store/slice/InputNode/InputNodeType'
import { FILE_TREE_TYPE_SET } from 'store/slice/FilesTree/FilesTreeType'
import { Box, TextField } from '@mui/material'
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
  }
  const filePath = useSelector(selectImageItemFilePath(itemId))

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
        onSelectFile={(path) => !Array.isArray(path) && onSelectImageFile(path)}
        onUploadFile={onUploadFileHandle}
        fileTreeType={FILE_TREE_TYPE_SET.IMAGE}
        selectButtonLabel="Select Image"
      />
      <StartEndIndex />
      <Showticklabels />
      <ShowLine />
      <ShowGrid />
      <ShowScale />
      <Zsmooth />
      <GradientColorPicker colors={colors} dispatchSetColor={dispathSetColor} />
      <div>
        <h3>Roi Setting</h3>
        <FilePathSelect
          selectedFilePath={roiItemFilePath}
          selectedNodeId={roiItemNodeId}
          onSelect={onSelectRoiFilePath}
          dataType={DATA_TYPE_SET.ROI}
          label={'Select Roi'}
        />
        <RoiAlpha />
      </div>
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
  const handleChange = (event: SelectChangeEvent<string | boolean>) => {
    dispatch(setImageItemZsmooth({ itemId, zsmooth: event.target.value }))
  }
  return (
    <FormControl variant="standard" sx={{ m: 1, minWidth: 120 }}>
      <InputLabel>smooth</InputLabel>
      <Select label="smooth" value={zsmooth} onChange={handleChange}>
        <MenuItem value={'best'}>best</MenuItem>
        <MenuItem value={'fast'}>fast</MenuItem>
        <MenuItem value={'false'}>False</MenuItem>
      </Select>
    </FormControl>
  )
}

const RoiAlpha: React.FC = () => {
  const itemId = React.useContext(SelectedItemIdContext)
  const dispatch = useDispatch()
  const roiAlpha = useSelector(selectImageItemRoiAlpha(itemId))
  const inputError = !(roiAlpha > 0)
  const onChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    const newValue = event.target.value === '' ? '' : Number(event.target.value)
    if (typeof newValue === 'number') {
      dispatch(setImageItemRoiAlpha({ itemId, roiAlpha: newValue }))
    }
  }
  return (
    <>
      <TextField
        error={inputError}
        type="number"
        inputProps={{
          step: 0.1,
          min: 0,
          max: 1.0,
        }}
        InputLabelProps={{
          shrink: true,
        }}
        onChange={onChange}
        value={roiAlpha}
        helperText={inputError ? 'index > 0' : undefined}
      />
      alpha
    </>
  )
}

const StartEndIndex: React.FC = () => {
  const itemId = React.useContext(SelectedItemIdContext)
  const [startIndex, onChangeStartIndex] = React.useState(
    useSelector(selectImageItemStartIndex(itemId)),
  )
  const [endIndex, onChangeEndIndex] = React.useState(
    useSelector(selectImageItemEndIndex(itemId)),
  )
  const inputError = !(startIndex > 0)
  const dispatch = useDispatch()
  const onStartChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    const newValue = event.target.value === '' ? '' : Number(event.target.value)
    if (typeof newValue === 'number') {
      onChangeStartIndex(newValue)
    }
  }
  const onEndChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    const newValue = event.target.value === '' ? '' : Number(event.target.value)
    if (typeof newValue === 'number') {
      onChangeEndIndex(newValue)
    }
  }

  const filePath = useSelector(selectImageItemFilePath(itemId))
  const onClickButton = () => {
    if (startIndex > 0) {
      dispatch(setImageItemStartIndex({ itemId, startIndex }))
      dispatch(setImageItemEndIndex({ itemId, endIndex }))
      dispatch(resetImageActiveIndex({ itemId, startIndex, endIndex }))
      if (filePath !== null) {
        dispatch(
          getImageData({
            path: filePath,
            startIndex: startIndex ?? 1,
            endIndex: endIndex ?? 10,
          }),
        )
      }
    }
  }

  return (
    <Box sx={{ display: 'flex', alignItems: 'flex-start' }}>
      <TextField
        error={inputError}
        type="number"
        inputProps={{
          step: 1,
          min: 0,
        }}
        InputLabelProps={{
          shrink: true,
        }}
        onChange={onStartChange}
        value={startIndex}
        helperText={inputError ? 'index > 0' : undefined}
      />
      ~
      <TextField
        type="number"
        InputLabelProps={{
          shrink: true,
        }}
        onChange={onEndChange}
        value={endIndex}
      />
      <Button
        size="small"
        className="ctrl_btn"
        variant="contained"
        onClick={onClickButton}
      >
        load
      </Button>
    </Box>
  )
}
