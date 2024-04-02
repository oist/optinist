import { memo, MouseEvent, useState } from "react"

import DoneIcon from "@mui/icons-material/Done"
import ErrorOutlineIcon from "@mui/icons-material/ErrorOutline"
import HorizontalRuleIcon from "@mui/icons-material/HorizontalRule"
import { IconButton, Popover, Typography } from "@mui/material"

import { EXPERIMENTS_STATUS } from "store/slice/Experiments/ExperimentsType"

interface ExperimentStatusIconProps {
  status: EXPERIMENTS_STATUS
  message?: string
}

export const ExperimentStatusIcon = memo(function ExperimentStatusIcon({
  status,
  message,
}: ExperimentStatusIconProps) {
  switch (status) {
    case "error":
      return <ErrorIcon message={message} />
    case "success":
      return <DoneIcon color="success" />
    case "running":
      return <HorizontalRuleIcon color="inherit" />
  }
})

interface ErrorIconProps {
  message?: string
}

const ErrorIcon = memo(function ErrorIcon({ message }: ErrorIconProps) {
  const [anchorEl, setAnchorEl] = useState<HTMLButtonElement | null>(null)
  const handleClick = (event: MouseEvent<HTMLButtonElement>) => {
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
          whiteSpace="pre-wrap"
        >
          {message}
        </Typography>
      </Popover>
    </>
  )
})
