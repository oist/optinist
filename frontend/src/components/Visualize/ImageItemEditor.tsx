import React from 'react'

import { useDispatch, useSelector } from 'react-redux'
import Switch from '@material-ui/core/Switch'
import FormControlLabel from '@material-ui/core/FormControlLabel'

import { selectImageItemShowticklabels } from 'store/slice/VisualizeItem/VisualizeItemSelectors'
import { SelectedItemIdContext } from './VisualizeItemEditor'
import { setImageItemShowticklabels } from 'store/slice/VisualizeItem/VisualizeItemSlice'

export const ImageItemEditor: React.FC = () => {
  return (
    <div>
      <Showticklabels />
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
