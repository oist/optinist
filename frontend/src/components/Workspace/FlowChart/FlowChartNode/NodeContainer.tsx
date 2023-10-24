import { FC, ReactNode } from "react"

import { Box, Typography } from "@mui/material"
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
      <Typography
        marginX={1}
        marginTop="auto"
        marginBottom="4px"
        color={grey[600]}
        fontSize={13}
        paragraph
        variant="body2"
        sx={{ overflowWrap: "break-word" }}
      >
        {nodeId}
      </Typography>
    </div>
  )
}
