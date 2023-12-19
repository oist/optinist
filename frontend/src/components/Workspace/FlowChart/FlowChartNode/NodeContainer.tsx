import { FC, ReactNode } from "react"

import { Box, Tooltip, Typography } from "@mui/material"
import { grey } from "@mui/material/colors"
import { alpha, useTheme } from "@mui/material/styles"

export const NodeContainer: FC<{
  children: ReactNode
  nodeId: string
  selected: boolean
  updated?: boolean
}> = ({ children, nodeId, selected, updated }) => {
  const theme = useTheme()

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
      <Tooltip title={nodeId} placement="bottom-start">
        <Typography
          marginX={1}
          marginTop="auto"
          marginBottom="4px"
          color={grey[600]}
          fontSize={14}
          variant="body2"
          overflow="ellipsis"
        >
          {nodeId}
        </Typography>
      </Tooltip>
    </div>
  )
}
