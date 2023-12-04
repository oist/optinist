import { ChangeEvent, FC, useContext } from "react"
import { useDispatch, useSelector } from "react-redux"

import MenuItem from "@mui/material/MenuItem"
import { SelectChangeEvent } from "@mui/material/Select"

import { ParamSection } from "components/common/ParamSection"
import { ParamSelect } from "components/common/ParamSelect"
import { ParamTextField } from "components/common/ParamTextField"
import { SelectedItemIdContext } from "components/Workspace/Visualize/VisualizeItemEditor"
import {
  selectVisualizeSaveFilename,
  selectVisualizeSaveFormat,
} from "store/slice/VisualizeItem/VisualizeItemSelectors"
import {
  setSaveFileName,
  setSaveFormat,
} from "store/slice/VisualizeItem/VisualizeItemSlice"

import "react-linear-gradient-picker/dist/index.css"

export const SaveFig: FC = () => {
  const itemId = useContext(SelectedItemIdContext)
  const saveFileName = useSelector(selectVisualizeSaveFilename(itemId))
  const saveFormat = useSelector(selectVisualizeSaveFormat(itemId))
  const dispatch = useDispatch()
  const handleChange = (event: SelectChangeEvent<string>) => {
    dispatch(setSaveFormat({ itemId, saveFormat: event.target.value }))
  }
  const onChangeFileName = (event: ChangeEvent<HTMLInputElement>) => {
    dispatch(setSaveFileName({ itemId, saveFileName: event.target.value }))
  }

  return (
    <ParamSection title="SaveFig">
      <ParamSelect
        label={"File Format"}
        value={saveFormat}
        onChange={handleChange}
      >
        <MenuItem value={"svg"}>svg</MenuItem>
        <MenuItem value={"png"}>png</MenuItem>
        <MenuItem value={"jpeg"}>jpeg</MenuItem>
        <MenuItem value={"webp"}>webp</MenuItem>
      </ParamSelect>
      <ParamTextField
        type="text"
        label={"Filename"}
        value={saveFileName}
        onChange={onChangeFileName}
      />
    </ParamSection>
  )
}
