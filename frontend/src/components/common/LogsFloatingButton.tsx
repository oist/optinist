import { memo } from "react"
import { useDispatch } from "react-redux"

import TerminalIcon from "@mui/icons-material/Terminal"
import { Box, Fab, Tooltip } from "@mui/material"
import { styled } from "@mui/material/styles"

import { openLogsModal } from "store/slice/LogsModal/LogsModalSlice"
import { AppDispatch } from "store/store"

export const LogsFloatingButton = memo(function LogsFloatingButton() {
  const dispatch = useDispatch<AppDispatch>()

  const handleClick = () => {
    dispatch(openLogsModal())
  }

  return (
    <FloatingButtonContainer>
      <Tooltip title="Show Logs" placement="right">
        <StyledFab onClick={handleClick} size="medium" aria-label="show logs">
          <TerminalIcon style={{ fontSize: 22 }} />
        </StyledFab>
      </Tooltip>
    </FloatingButtonContainer>
  )
})

const FloatingButtonContainer = styled(Box)({
  position: "fixed",
  bottom: 20,
  left: 20,
  zIndex: 1000,
})

const StyledFab = styled(Fab)(({ theme }) => ({
  backgroundColor: "#ffffff",
  color: theme.palette.primary.main,
  boxShadow: theme.shadows[4],
  border: `1px solid ${theme.palette.divider}`,
  "&:hover": {
    backgroundColor: "#f5f5f5",
    boxShadow: theme.shadows[8],
  },
}))
