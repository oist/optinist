import { FC, useContext } from "react"
import { useSelector, useDispatch } from "react-redux"

import { Grid } from "@mui/material"

import { FieldLabel, ParamSection } from "components/common/ParamSection"
import { ParamSwitch } from "components/common/ParamSwitch"
import { GradientColorPicker } from "components/Workspace/Visualize/Editor/GradientColorPicker"
import { SaveFig } from "components/Workspace/Visualize/Editor/SaveFig"
import { SelectedItemIdContext } from "components/Workspace/Visualize/VisualizeItemEditor"
import {
  selectHeatMapItemColors,
  selectHeatMapItemShowScale,
} from "store/slice/VisualizeItem/VisualizeItemSelectors"
import {
  setHeatMapItemColors,
  setHeatMapItemShowScale,
} from "store/slice/VisualizeItem/VisualizeItemSlice"
import { ColorType } from "store/slice/VisualizeItem/VisualizeItemType"

export const HeatmapItemEditor: FC = () => {
  const itemId = useContext(SelectedItemIdContext)
  const dispatch = useDispatch()
  const colors = useSelector(selectHeatMapItemColors(itemId))
  const dispathSetColor = (colorCode: ColorType[]) => {
    dispatch(setHeatMapItemColors({ itemId, colors: colorCode }))
  }
  return (
    <>
      <ParamSection title="Heatmap">
        <ShowScale />
        <Grid container component="label" alignItems="center">
          <Grid item xs={8}>
            <FieldLabel>Pick Color</FieldLabel>
          </Grid>
          <Grid item xs={4}>
            <GradientColorPicker
              colors={colors}
              dispatchSetColor={dispathSetColor}
            />
          </Grid>
        </Grid>
      </ParamSection>
      <SaveFig />
    </>
  )
}

const ShowScale: FC = () => {
  const itemId = useContext(SelectedItemIdContext)
  const showscale = useSelector(selectHeatMapItemShowScale(itemId))
  const dispatch = useDispatch()
  const toggleChecked = () => {
    dispatch(setHeatMapItemShowScale({ itemId, showscale: !showscale }))
  }
  return (
    <ParamSwitch label="ShowScale" value={showscale} onChange={toggleChecked} />
  )
}
