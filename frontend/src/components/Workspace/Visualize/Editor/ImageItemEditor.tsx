import { ChangeEvent, FC, useContext, useState } from "react"
import { useDispatch, useSelector } from "react-redux"

import "react-linear-gradient-picker/dist/index.css"
import { Box, Grid, TextField, Typography } from "@mui/material"
import Button from "@mui/material/Button"
import MenuItem from "@mui/material/MenuItem"
import { SelectChangeEvent } from "@mui/material/Select"

import { ParamSection, FieldLabel } from "components/common/ParamSection"
import { ParamSelect } from "components/common/ParamSelect"
import { ParamSwitch } from "components/common/ParamSwitch"
import { ParamTextField } from "components/common/ParamTextField"
import { GradientColorPicker } from "components/Workspace/Visualize/Editor/GradientColorPicker"
import { SaveFig } from "components/Workspace/Visualize/Editor/SaveFig"
import { SelectedItemIdContext } from "components/Workspace/Visualize/VisualizeItemEditor"
import { getImageData } from "store/slice/DisplayData/DisplayDataActions"
import { getLabelByPath } from "store/slice/FlowElement/FlowElementUtils"
import {
  selectImageItemShowGrid,
  selectImageItemShowLine,
  selectImageItemShowticklabels,
  selectImageItemZsmooth,
  selectImageItemShowScale,
  selectImageItemColors,
  selectImageItemStartIndex,
  selectImageItemEndIndex,
  selectImageItemRoiAlpha,
  selectImageItemFilePath,
  selectImageItemAlpha,
} from "store/slice/VisualizeItem/VisualizeItemSelectors"
import {
  setImageItemShowGrid,
  setImageItemShowLine,
  setImageItemShowticklabels,
  setImageItemZsmooth,
  setImageItemShowScale,
  setImageItemStartIndex,
  setImageItemEndIndex,
  setImageItemColors,
  resetImageActiveIndex,
  setImageItemRoiAlpha,
  setImageItemAlpha,
} from "store/slice/VisualizeItem/VisualizeItemSlice"
import { ColorType } from "store/slice/VisualizeItem/VisualizeItemType"
import { selectCurrentWorkspaceId } from "store/slice/Workspace/WorkspaceSelector"
import { AppDispatch } from "store/store"

export const ImageItemEditor: FC = () => {
  const itemId = useContext(SelectedItemIdContext)
  const dispatch = useDispatch()
  const filePath = useSelector(selectImageItemFilePath(itemId))
  const colors = useSelector(selectImageItemColors(itemId))
  const dispathSetColor = (colorCode: ColorType[]) => {
    dispatch(setImageItemColors({ itemId, colors: colorCode }))
  }

  return (
    <>
      <ParamSection title="Image">
        <Typography marginBottom={2}>
          {filePath ? getLabelByPath(filePath) : "No file is selected"}
        </Typography>
        <StartEndIndex />
        <Showticklabels />
        <ShowLine />
        <ShowGrid />
        <ShowScale />
        <Zsmooth />
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
        <Alpha />
        <RoiAlpha />
      </ParamSection>
      <SaveFig />
    </>
  )
}

const Showticklabels: FC = () => {
  const itemId = useContext(SelectedItemIdContext)
  const showticklabels = useSelector(selectImageItemShowticklabels(itemId))
  const dispatch = useDispatch()
  const toggleChecked = () => {
    dispatch(
      setImageItemShowticklabels({ itemId, showticklabels: !showticklabels }),
    )
  }
  return (
    <ParamSwitch
      label={"ShowTickLabels"}
      value={showticklabels}
      onChange={toggleChecked}
    />
  )
}

const ShowLine: FC = () => {
  const itemId = useContext(SelectedItemIdContext)
  const showline = useSelector(selectImageItemShowLine(itemId))
  const dispatch = useDispatch()
  const toggleChecked = () => {
    dispatch(setImageItemShowLine({ itemId, showline: !showline }))
  }
  return (
    <ParamSwitch label={"ShowLine"} value={showline} onChange={toggleChecked} />
  )
}

const ShowGrid: FC = () => {
  const itemId = useContext(SelectedItemIdContext)
  const showgrid = useSelector(selectImageItemShowGrid(itemId))
  const dispatch = useDispatch()
  const toggleChecked = () => {
    dispatch(setImageItemShowGrid({ itemId, showgrid: !showgrid }))
  }
  return (
    <ParamSwitch label={"ShowGrid"} value={showgrid} onChange={toggleChecked} />
  )
}

const ShowScale: FC = () => {
  const itemId = useContext(SelectedItemIdContext)
  const showscale = useSelector(selectImageItemShowScale(itemId))
  const dispatch = useDispatch()
  const toggleChecked = () => {
    dispatch(setImageItemShowScale({ itemId, showscale: !showscale }))
  }
  return (
    <ParamSwitch
      label={"ShowScale"}
      value={showscale}
      onChange={toggleChecked}
    />
  )
}

const Zsmooth: FC = () => {
  const itemId = useContext(SelectedItemIdContext)
  const zsmooth = useSelector(selectImageItemZsmooth(itemId))
  const dispatch = useDispatch()
  const handleChange = (event: SelectChangeEvent<string | boolean>) => {
    dispatch(setImageItemZsmooth({ itemId, zsmooth: event.target.value }))
  }
  return (
    <ParamSelect
      label="Smooth"
      onChange={handleChange}
      value={zsmooth as string}
    >
      <MenuItem value={"best"}>best</MenuItem>
      <MenuItem value={"fast"}>fast</MenuItem>
      <MenuItem value={"false"}>False</MenuItem>
    </ParamSelect>
  )
}

const Alpha: FC = () => {
  const itemId = useContext(SelectedItemIdContext)
  const dispatch = useDispatch()
  const alpha = useSelector(selectImageItemAlpha(itemId))
  const inputError = !(alpha > 0)
  const onChange = (event: ChangeEvent<HTMLInputElement>) => {
    const newValue = event.target.value === "" ? "" : Number(event.target.value)
    if (typeof newValue === "number") {
      dispatch(setImageItemAlpha({ itemId, alpha: newValue }))
    }
  }
  return (
    <ParamTextField
      label={"Image Alpha"}
      type="number"
      value={alpha}
      inputProps={{
        step: 0.1,
        min: 0,
        max: 1.0,
      }}
      onChange={onChange}
      error={inputError}
      helperText={inputError ? "index > 0" : undefined}
    />
  )
}

const RoiAlpha: FC = () => {
  const itemId = useContext(SelectedItemIdContext)
  const dispatch = useDispatch()
  const roiAlpha = useSelector(selectImageItemRoiAlpha(itemId))
  const inputError = !(roiAlpha > 0)
  const onChange = (event: ChangeEvent<HTMLInputElement>) => {
    const newValue = event.target.value === "" ? "" : Number(event.target.value)
    if (typeof newValue === "number") {
      dispatch(setImageItemRoiAlpha({ itemId, roiAlpha: newValue }))
    }
  }
  return (
    <ParamTextField
      label={"Roi Alpha"}
      type="number"
      value={roiAlpha}
      inputProps={{
        step: 0.1,
        min: 0,
        max: 1.0,
      }}
      onChange={onChange}
      error={inputError}
      helperText={inputError ? "index > 0" : undefined}
    />
  )
}

const StartEndIndex: FC = () => {
  const workspaceId = useSelector(selectCurrentWorkspaceId)
  const itemId = useContext(SelectedItemIdContext)
  const [startIndex, onChangeStartIndex] = useState(
    useSelector(selectImageItemStartIndex(itemId)),
  )
  const [endIndex, onChangeEndIndex] = useState(
    useSelector(selectImageItemEndIndex(itemId)),
  )
  const inputError = !(startIndex > 0)
  const dispatch = useDispatch<AppDispatch>()
  const onStartChange = (event: ChangeEvent<HTMLInputElement>) => {
    const newValue = event.target.value === "" ? "" : Number(event.target.value)
    if (typeof newValue === "number") {
      onChangeStartIndex(newValue)
    }
  }
  const onEndChange = (event: ChangeEvent<HTMLInputElement>) => {
    const newValue = event.target.value === "" ? "" : Number(event.target.value)
    if (typeof newValue === "number") {
      onChangeEndIndex(newValue)
    }
  }

  const filePath = useSelector(selectImageItemFilePath(itemId))
  const onClickButton = () => {
    if (startIndex > 0) {
      dispatch(setImageItemStartIndex({ itemId, startIndex }))
      dispatch(setImageItemEndIndex({ itemId, endIndex }))
      dispatch(resetImageActiveIndex({ itemId, startIndex, endIndex }))
      if (workspaceId && filePath !== null) {
        dispatch(
          getImageData({
            workspaceId,
            path: filePath,
            startIndex: startIndex ?? 1,
            endIndex: endIndex ?? 10,
          }),
        )
      }
    }
  }

  return (
    <Box marginBottom={2}>
      <FieldLabel>Start/End Index</FieldLabel>
      <Box sx={{ display: "flex", alignItems: "flex-start" }}>
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
          helperText={inputError ? "index > 0" : undefined}
          style={{ marginRight: 8 }}
        />
        ~
        <TextField
          type="number"
          InputLabelProps={{
            shrink: true,
          }}
          onChange={onEndChange}
          value={endIndex}
          style={{ marginLeft: 8, marginRight: 8 }}
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
    </Box>
  )
}
