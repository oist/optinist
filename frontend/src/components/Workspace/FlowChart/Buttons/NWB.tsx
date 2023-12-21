import { memo } from "react"
import { useDispatch, useSelector } from "react-redux"

import TuneIcon from "@mui/icons-material/Tune"
import { IconButton, Tooltip } from "@mui/material"

import { selectPipelineIsStartedSuccess } from "store/slice/Pipeline/PipelineSelectors"
import { toggleNwb } from "store/slice/RightDrawer/RightDrawerSlice"

export const NWBSettingButton = memo(function NWBSettingButton() {
  const dispatch = useDispatch()
  const handleClick = () => {
    dispatch(toggleNwb())
  }
  const isPending = useSelector(selectPipelineIsStartedSuccess)
  return (
    <Tooltip title="NWB settings">
      <IconButton onClick={handleClick} color="primary" disabled={!!isPending}>
        <TuneIcon />
      </IconButton>
    </Tooltip>
  )
})
