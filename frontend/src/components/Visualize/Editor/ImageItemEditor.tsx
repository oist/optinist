import React, { useState } from 'react'

import { useDispatch, useSelector } from 'react-redux'
import Switch from '@material-ui/core/Switch'
import FormControlLabel from '@material-ui/core/FormControlLabel'
import Select from '@material-ui/core/Select'
import MenuItem from '@material-ui/core/MenuItem'

import {
  selectImageItemShowGrid,
  selectImageItemShowLine,
  selectImageItemShowticklabels,
  selectImageItemZsmooth,
  selectImageItemShowScale,
  selectImageItemColors,
  selectImageItemMaxIndex,
  selectVisualizeDataFilePath,
} from 'store/slice/VisualizeItem/VisualizeItemSelectors'
import { SelectedItemIdContext } from '../VisualizeItemEditor'

import {
  setImageItemShowGrid,
  setImageItemShowLine,
  setImageItemShowticklabels,
  setImageItemZsmooth,
  setImageItemShowScale,
  setImageItemColors,
  setDisplayDataPath,
  setImageItemMaxIndex,
} from 'store/slice/VisualizeItem/VisualizeItemSlice'

import {
  GradientPicker,
  GradientPickerPopover,
} from 'react-linear-gradient-picker'
import { SketchPicker } from 'react-color'
import { PALETTE_COLOR_SHAPE_TYPE } from 'react-linear-gradient-picker'
import 'react-linear-gradient-picker/dist/index.css'
import { FileSelect } from 'components/FlowChart/FlowChartNode/FileSelect'
import { useFileUploader } from 'store/slice/FileUploader/FileUploaderHook'
import { FILE_TYPE_SET } from 'store/slice/InputNode/InputNodeType'
import { FILE_TREE_TYPE_SET } from 'store/slice/FilesTree/FilesTreeType'
import { TextField } from '@material-ui/core'

export const ImageItemEditor: React.FC = () => {
  const itemId = React.useContext(SelectedItemIdContext)
  // const [filePath, setFilePath] = useState('')
  const dispatch = useDispatch()
  const onSelectFile = (path: string) => {
    // setFilePath(path)
    dispatch(setDisplayDataPath({ nodeId: null, filePath: path, itemId }))
  }
  const filePath = useSelector(selectVisualizeDataFilePath(itemId))

  const {
    // filePath: uploadedFilePath,
    onUploadFile,
    // pending,
    // uninitialized,
    // progress,
    // error,
  } = useFileUploader(FILE_TYPE_SET.IMAGE)
  const onUploadFileHandle = (formData: FormData, fileName: string) => {
    onUploadFile(formData, fileName)
  }

  return (
    <div style={{ margin: '10px' }}>
      <FileSelect
        filePath={filePath ?? ''}
        onSelectFile={onSelectFile}
        onUploadFile={onUploadFileHandle}
        fileTreeType={FILE_TREE_TYPE_SET.IMAGE}
        selectButtonLabel="Select Image"
      />
      <MaxIndex />
      <Showticklabels />
      <ShowLine />
      <ShowGrid />
      <ShowScale />
      <Zsmooth />
      <GradientColorPicker />
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

const MaxIndex: React.FC = () => {
  const itemId = React.useContext(SelectedItemIdContext)
  const maxIndex = useSelector(selectImageItemMaxIndex(itemId))
  const dispatch = useDispatch()
  const onChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    const newValue = event.target.value === '' ? '' : Number(event.target.value)
    if (typeof newValue === 'number') {
      dispatch(setImageItemMaxIndex({ itemId, maxIndex: newValue }))
    }
  }
  return (
    <FormControlLabel
      control={
        <TextField
          type="number"
          InputLabelProps={{
            shrink: true,
          }}
          onChange={onChange}
          defaultValue={maxIndex}
        />
      }
      label="maxIndex"
    />
  )
}

const GradientColorPicker: React.FC = () => {
  const itemId = React.useContext(SelectedItemIdContext)
  const colors = useSelector(selectImageItemColors(itemId))
  const dispatch = useDispatch()

  const palette: PALETTE_COLOR_SHAPE_TYPE[] = colors.map((value) => {
    return {
      offset: value.offset,
      color: value.rgb,
    }
  })

  const onPaletteChange = (palette: PALETTE_COLOR_SHAPE_TYPE[]) => {
    const colorCode = palette.map((value) => {
      const color = value.color
      const rgbStr = color.replace(/[^0-9,]/g, '').split(',')
      const rgb = {
        r: Number(rgbStr[0]),
        g: Number(rgbStr[1]),
        b: Number(rgbStr[2]),
      }
      return {
        rgb: `rgb(${rgb.r}, ${rgb.g}, ${rgb.b})`,
        offset: value.offset,
      }
    })
    dispatch(setImageItemColors({ itemId, colors: colorCode }))
  }

  const [open, setOpen] = useState(false)

  return (
    // <GradientPicker
    //   width={200}
    //   maxStops={10}
    //   paletteHeight={32}
    //   palette={palette}
    //   onPaletteChange={onPaletteChange}
    //   flatStyle={true}
    // >
    //   <WrappedSketchPicker />
    // </GradientPicker>
    <GradientPickerPopover
      open={open}
      setOpen={() => setOpen(!open)}
      // showAnglePicker={true}
      width={150}
      maxStops={10}
      paletteHeight={25}
      palette={palette}
      onPaletteChange={onPaletteChange}
      flatStyle={true}
    >
      <WrappedSketchPicker />
    </GradientPickerPopover>
  )
}

type WrapperPropTypes = {
  onSelect?: (color: string, opacity?: number) => void
  color?: string
}

const WrappedSketchPicker: React.FC<WrapperPropTypes> = ({
  onSelect,
  color,
}) => {
  return (
    <SketchPicker
      color={color}
      width="150px"
      // styles={{width: "10px"}}
      onChange={(c) => {
        const { r, g, b, a } = c.rgb
        onSelect?.(`rgb(${r}, ${g}, ${b})`, a)
      }}
    />
  )
}
