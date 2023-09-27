import React from 'react'
import { useDispatch, useSelector } from 'react-redux'
import Select, { SelectChangeEvent } from '@mui/material/Select'
import InputLabel from '@mui/material/InputLabel'
import MenuItem from '@mui/material/MenuItem'
import FormControl from '@mui/material/FormControl'
import TextField from '@mui/material/TextField'

import {
  selectVisualizeSaveFilename,
  selectVisualizeSaveFormat,
} from 'store/slice/VisualizeItem/VisualizeItemSelectors'
import { SelectedItemIdContext } from '../VisualizeItemEditor'

import {
  setSaveFileName,
  setSaveFormat,
} from 'store/slice/VisualizeItem/VisualizeItemSlice'

import 'react-linear-gradient-picker/dist/index.css'

export const SaveFig: React.FC = () => {
  const itemId = React.useContext(SelectedItemIdContext)
  const saveFileName = useSelector(selectVisualizeSaveFilename(itemId))
  const saveFormat = useSelector(selectVisualizeSaveFormat(itemId))
  const dispatch = useDispatch()
  const handleChange = (event: SelectChangeEvent<string>) => {
    dispatch(setSaveFormat({ itemId, saveFormat: event.target.value }))
  }
  const onChangeFileName = (event: React.ChangeEvent<HTMLInputElement>) => {
    dispatch(setSaveFileName({ itemId, saveFileName: event.target.value }))
  }

  return (
    <>
      <h3>SaveFig</h3>
      <FormControl
        variant="standard"
        sx={{ minWidth: 120, width: '100%', marginBottom: 1 }}
      >
        <InputLabel>format</InputLabel>
        <Select label="smooth" value={saveFormat} onChange={handleChange}>
          <MenuItem value={'svg'}>svg</MenuItem>
          <MenuItem value={'png'}>png</MenuItem>
          <MenuItem value={'jpeg'}>jpeg</MenuItem>
          <MenuItem value={'webp'}>webp</MenuItem>
        </Select>
      </FormControl>
      <TextField
        style={{ width: '100%' }}
        label={'Filename'}
        InputLabelProps={{
          shrink: true,
        }}
        onChange={onChangeFileName}
        value={saveFileName}
      />
    </>
  )
}
