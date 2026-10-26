import { ChangeEvent, FC, useContext } from "react"
import { useSelector, useDispatch } from "react-redux"

import { Typography } from "@mui/material"

import { FileNameChip } from "components/common/ParamSection"
import { ParamSwitch } from "components/common/ParamSwitch"
import { ParamTextField } from "components/common/ParamTextField"
import { SelectedItemIdContext } from "components/Workspace/Visualize/VisualizeItemEditor"
import {
  selectCsvItemSetHeader,
  selectCsvItemSetIndex,
  selectCsvItemTranspose,
  selectVisualizeDataFilePath,
} from "store/slice/VisualizeItem/VisualizeItemSelectors"
import {
  setCsvItemSetHeader,
  setCsvItemSetIndex,
  setCsvItemTranspose,
} from "store/slice/VisualizeItem/VisualizeItemSlice"

export const CsvItemEditor: FC = () => {
  const itemId = useContext(SelectedItemIdContext)
  const filePath = useSelector(selectVisualizeDataFilePath(itemId))

  return (
    <>
      <Typography variant="h6" fontWeight="bold">
        Csv
      </Typography>
      <FileNameChip filePath={filePath} />
      <Transpose />
      <SetHeader />
      <SetIndex />
    </>
  )
}

const Transpose: FC = () => {
  const itemId = useContext(SelectedItemIdContext)
  const transpose = useSelector(selectCsvItemTranspose(itemId))
  const dispatch = useDispatch()
  const toggleChecked = () => {
    dispatch(setCsvItemTranspose({ itemId, transpose: !transpose }))
  }
  return (
    <ParamSwitch
      label={"Transpose"}
      value={transpose}
      onChange={toggleChecked}
    />
  )
}

const SetHeader: FC = () => {
  const itemId = useContext(SelectedItemIdContext)
  const setHeader = useSelector(selectCsvItemSetHeader(itemId))

  const dispatch = useDispatch()
  const onChangeSetHeader = (event: ChangeEvent<HTMLInputElement>) => {
    const newValue =
      event.target.value === "" ? null : Number(event.target.value)
    if (newValue === null || newValue >= 0) {
      dispatch(setCsvItemSetHeader({ itemId, setHeader: newValue }))
    }
  }

  return (
    <ParamTextField
      label="Header"
      value={setHeader ?? ""}
      type="number"
      onChange={onChangeSetHeader}
    />
  )
}

const SetIndex: FC = () => {
  const itemId = useContext(SelectedItemIdContext)
  const setIndex = useSelector(selectCsvItemSetIndex(itemId))
  const dispatch = useDispatch()
  const toggleChecked = () => {
    dispatch(setCsvItemSetIndex({ itemId, setIndex: !setIndex }))
  }
  return (
    <ParamSwitch label={"SetIndex"} value={setIndex} onChange={toggleChecked} />
  )
}
