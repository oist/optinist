import React from 'react'
import { useSelector, useDispatch } from 'react-redux'
import { FormControlLabel, Switch } from '@material-ui/core'
import {
  selectHeatMapItemColors,
  selectHeatMapItemShowScale,
} from 'store/slice/VisualizeItem/VisualizeItemSelectors'
import {
  setHeatMapItemColors,
  setHeatMapItemShowScale,
} from 'store/slice/VisualizeItem/VisualizeItemSlice'
import { SelectedItemIdContext } from '../VisualizeItemEditor'
import { GradientColorPicker } from './GradientColorPicker'
import { ColorType } from 'store/slice/VisualizeItem/VisualizeItemType'

export const HeatmapItemEditor: React.FC = () => {
  const itemId = React.useContext(SelectedItemIdContext)
  const dispatch = useDispatch()
  const colors = useSelector(selectHeatMapItemColors(itemId))
  const dispathSetColor = (colorCode: ColorType[]) => {
    dispatch(setHeatMapItemColors({ itemId, colors: colorCode }))
  }
  return (
    <>
      <ShowScale />
      <GradientColorPicker colors={colors} dispatchSetColor={dispathSetColor} />
    </>
  )
}

const ShowScale: React.FC = () => {
  const itemId = React.useContext(SelectedItemIdContext)
  const showscale = useSelector(selectHeatMapItemShowScale(itemId))
  const dispatch = useDispatch()
  const toggleChecked = () => {
    dispatch(setHeatMapItemShowScale({ itemId, showscale: !showscale }))
  }
  return (
    <FormControlLabel
      control={<Switch checked={showscale} onChange={toggleChecked} />}
      label="showscale"
    />
  )
}
