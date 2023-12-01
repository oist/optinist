import { FC, ReactNode } from "react"

import { useSnackbar } from "notistack"

import ContentCopyRoundedIcon from "@mui/icons-material/ContentCopyRounded"
import { Box, IconButton, Tooltip, Typography } from "@mui/material"
import { grey } from "@mui/material/colors"
import { alpha, useTheme } from "@mui/material/styles"

export const NodeContainer: FC<{
  children: ReactNode
  nodeId: string
  selected: boolean
  updated?: boolean
}> = ({ children, nodeId, selected, updated }) => {
  const theme = useTheme()
  const { enqueueSnackbar } = useSnackbar()

  const handleNodeIdClick = () => {
    navigator.clipboard.writeText(nodeId)
    enqueueSnackbar("Node ID copied to clipboard", { variant: "success" })
  }

  let backgroundColor
  if (updated) {
    if (selected) {
      backgroundColor = alpha(theme.palette.warning.light, 0.3)
    } else {
      backgroundColor = alpha(theme.palette.warning.light, 0.1)
    }
  } else {
    if (selected) {
      backgroundColor = alpha(theme.palette.primary.light, 0.2)
    }
  }

  return (
    <div
      style={{
        height: "100%",
        width: "100%",
        background: backgroundColor,
        display: "flex",
        flexDirection: "column",
      }}
    >
      <Box margin={1}>{children}</Box>
      <div
        style={{
          marginLeft: 8,
          marginRight: 8,
          marginTop: "auto",
          display: "flex",
          flexDirection: "row",
          alignItems: "center",
        }}
      >
        <Typography
          color={grey[600]}
          fontSize={13}
          paragraph
          variant="body2"
          marginBottom={0}
          sx={{ overflowWrap: "break-word" }}
        >
          {nodeId}
        </Typography>
        <Tooltip title="copy node ID to clipboard">
          <IconButton onClick={handleNodeIdClick} size="small">
            <ContentCopyRoundedIcon fontSize="small" />
          </IconButton>
        </Tooltip>
      </div>
    </div>
  )
}
