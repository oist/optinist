import { ChangeEvent, FC, useContext, useState } from "react"
import { useDispatch, useSelector } from "react-redux"

import "react-linear-gradient-picker/dist/index.css"
import { Box, TextField } from "@mui/material"
import Button from "@mui/material/Button"
import FormControl from "@mui/material/FormControl"
import FormControlLabel from "@mui/material/FormControlLabel"
import InputLabel from "@mui/material/InputLabel"
import MenuItem from "@mui/material/MenuItem"
import Select, { SelectChangeEvent } from "@mui/material/Select"
import Switch from "@mui/material/Switch"

import { FILE_TREE_TYPE_SET } from "api/files/Files"
import { FileSelectImple } from "components/Workspace/FlowChart/FlowChartNode/FileSelect"
import { GradientColorPicker } from "components/Workspace/Visualize/Editor/GradientColorPicker"
import { SaveFig } from "components/Workspace/Visualize/Editor/SaveFig"
import { SelectedItemIdContext } from "components/Workspace/Visualize/VisualizeItemEditor"
import { getImageData } from "store/slice/DisplayData/DisplayDataActions"
import { DATA_TYPE_SET } from "store/slice/DisplayData/DisplayDataType"
import { useFileUploader } from "store/slice/FileUploader/FileUploaderHook"
import { FILE_TYPE_SET } from "store/slice/InputNode/InputNodeType"
import { setNewDisplayDataPath } from "store/slice/VisualizeItem/VisualizeItemActions"
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
  selectDisplayDataIsSingle,
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

  const isSingleData = useSelector(selectDisplayDataIsSingle(itemId))
  const onSelectImageFile = (newPath: string) => {
    const basePayload = {
      itemId,
      nodeId: null,
      filePath: newPath,
    }
    dispatch(
      setNewDisplayDataPath(
        isSingleData && filePath != null
          ? {
              ...basePayload,
              deleteData: true,
              prevDataType: DATA_TYPE_SET.IMAGE,
              prevFilePath: filePath,
            }
          : {
              ...basePayload,
              deleteData: false,
            },
      ),
    )
  }

  const { onUploadFile } = useFileUploader({ fileType: FILE_TYPE_SET.IMAGE })
  const onUploadFileHandle = (formData: FormData, fileName: string) => {
    onUploadFile(formData, fileName)
  }

  const colors = useSelector(selectImageItemColors(itemId))
  const dispathSetColor = (colorCode: ColorType[]) => {
    dispatch(setImageItemColors({ itemId, colors: colorCode }))
  }

  return (
    <div style={{ margin: "10px", padding: 10 }}>
      <FileSelectImple
        filePath={filePath ?? ""}
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
      <Alpha />
      <RoiAlpha />
      <SaveFig />
    </div>
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
    <FormControlLabel
      control={<Switch checked={showticklabels} onChange={toggleChecked} />}
      label="Showticklabels"
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
    <FormControlLabel
      control={<Switch checked={showline} onChange={toggleChecked} />}
      label="ShowLine"
    />
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
    <FormControlLabel
      control={<Switch checked={showgrid} onChange={toggleChecked} />}
      label="ShowGrid"
    />
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
    <FormControlLabel
      control={<Switch checked={showscale} onChange={toggleChecked} />}
      label="ShowScale"
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
    <FormControl variant="standard" sx={{ m: 1, minWidth: 120 }}>
      <InputLabel>smooth</InputLabel>
      <Select label="smooth" value={zsmooth} onChange={handleChange}>
        <MenuItem value={"best"}>best</MenuItem>
        <MenuItem value={"fast"}>fast</MenuItem>
        <MenuItem value={"false"}>False</MenuItem>
      </Select>
    </FormControl>
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
    <>
      <TextField
        style={{ width: "100%" }}
        label={"image alpha"}
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
        value={alpha}
        helperText={inputError ? "index > 0" : undefined}
      />
    </>
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
    <>
      <TextField
        style={{ width: "100%" }}
        label={"roi alpha"}
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
        helperText={inputError ? "index > 0" : undefined}
      />
    </>
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
