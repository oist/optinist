import React from 'react'
import { useSelector, useDispatch } from 'react-redux'
import { FormControlLabel, TextField } from '@material-ui/core'
import {
  selectScatterItemXIndex,
  selectScatterItemYIndex,
} from 'store/slice/VisualizeItem/VisualizeItemSelectors'
import {
  setScatterItemXIndex,
  setScatterItemYIndex,
} from 'store/slice/VisualizeItem/VisualizeItemSlice'
import { SelectedItemIdContext } from '../VisualizeItemEditor'

export const ScatterItemEditor: React.FC = () => {
  return (
    <div style={{ margin: '10px', padding: 10 }}>
      <XIndex />
      <YIndex />
    </div>
  )
}

const XIndex: React.FC = () => {
  const itemId = React.useContext(SelectedItemIdContext)
  const xIndex = useSelector(selectScatterItemXIndex(itemId))

  const dispatch = useDispatch()
  const onChangeXIndex = (event: React.ChangeEvent<HTMLInputElement>) => {
    const newValue =
      event.target.value === '' ? null : Number(event.target.value)
    if (newValue !== null && newValue >= 0) {
      dispatch(setScatterItemXIndex({ itemId, xIndex: newValue }))
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
            onChange={onChangeXIndex}
            defaultValue={xIndex}
          />
          xIndex
        </>
      }
      label=""
    />
  )
}

const YIndex: React.FC = () => {
  const itemId = React.useContext(SelectedItemIdContext)
  const yIndex = useSelector(selectScatterItemYIndex(itemId))

  const dispatch = useDispatch()
  const onChangeYIndex = (event: React.ChangeEvent<HTMLInputElement>) => {
    const newValue =
      event.target.value === '' ? null : Number(event.target.value)
    if (newValue !== null && newValue >= 0) {
      dispatch(setScatterItemYIndex({ itemId, yIndex: newValue }))
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
            onChange={onChangeYIndex}
            defaultValue={yIndex}
          />
          yIndex
        </>
      }
      label=""
    />
  )
}
