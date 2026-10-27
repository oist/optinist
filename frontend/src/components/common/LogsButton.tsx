import { memo } from "react"
import { useDispatch } from "react-redux"

import TerminalIcon from "@mui/icons-material/Terminal"
import { IconButton, Tooltip } from "@mui/material"

import { openLogsModal } from "store/slice/LogsModal/LogsModalSlice"
import { AppDispatch } from "store/store"

export const LogsButton = memo(function LogsButton() {
  const dispatch = useDispatch<AppDispatch>()

  const handleClick = () => {
    dispatch(openLogsModal())
  }

  return (
    <Tooltip title="Show Logs">
      <IconButton onClick={handleClick} color={"primary"}>
        <TerminalIcon style={{ fontSize: 26 }} />
      </IconButton>
    </Tooltip>
  )
})
