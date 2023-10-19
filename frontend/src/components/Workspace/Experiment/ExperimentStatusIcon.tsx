import React from "react"

import DoneIcon from "@mui/icons-material/Done"
import ErrorOutlineIcon from "@mui/icons-material/ErrorOutline"
import HorizontalRuleIcon from "@mui/icons-material/HorizontalRule"
import { IconButton, Popover, Typography } from "@mui/material"

import { EXPERIMENTS_STATUS } from "store/slice/Experiments/ExperimentsType"

export const ExperimentStatusIcon = React.memo<{
  status: EXPERIMENTS_STATUS
  message?: string
}>(({ status, message }) => {
  switch (status) {
    case "error":
      return <ErrorIcon message={message} />
    case "success":
      return <DoneIcon color="success" />
    case "running":
      return <HorizontalRuleIcon color="inherit" />
  }
})

const ErrorIcon = React.memo<{ message?: string }>(({ message }) => {
  const [anchorEl, setAnchorEl] = React.useState<HTMLButtonElement | null>(null)
  const handleClick = (event: React.MouseEvent<HTMLButtonElement>) => {
    setAnchorEl(event.currentTarget)
  }
  const handleClose = () => {
    setAnchorEl(null)
  }

  const open = Boolean(anchorEl)
  const id = open ? "error-message-popover" : undefined

  return message == null ? (
    <ErrorOutlineIcon color="error" />
  ) : (
    <>
      <IconButton
        aria-describedby={id}
        onClick={handleClick}
        size="small"
        color="error"
        style={{ padding: 0 }}
      >
        <ErrorOutlineIcon color="error" />
      </IconButton>
      <Popover
        id={id}
        open={open}
        anchorEl={anchorEl}
        onClose={handleClose}
        anchorOrigin={{
          vertical: "bottom",
          horizontal: "center",
        }}
        slotProps={{ paper: { sx: { width: "60%" } } }}
      >
        <Typography
          variant="body2"
          paragraph
          color="error"
          padding={2}
          marginBottom={0}
        >
          {message}
        </Typography>
      </Popover>
    </>
  )
})
