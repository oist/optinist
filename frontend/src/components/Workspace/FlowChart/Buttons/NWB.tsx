import { memo } from "react"
import { useDispatch } from "react-redux"

import TuneIcon from "@mui/icons-material/Tune"
import { IconButton, Tooltip } from "@mui/material"

import { toggleNwb } from "store/slice/RightDrawer/RightDrawerSlice"

export const NWBSettingButton = memo(function NWBSettingButton() {
  const dispatch = useDispatch()
  const handleClick = () => {
    dispatch(toggleNwb())
  }
  return (
    <Tooltip title="NWB settings">
      <IconButton onClick={handleClick}>
        <TuneIcon color="primary" />
      </IconButton>
    </Tooltip>
  )
})
