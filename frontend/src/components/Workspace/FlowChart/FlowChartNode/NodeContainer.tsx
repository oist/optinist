import React from "react"
import { alpha, useTheme } from "@mui/material/styles"
import { Box, Typography } from "@mui/material"
import { grey } from "@mui/material/colors"

export const NodeContainer: React.FC<{
  children: React.ReactNode
  nodeId: string
  selected: boolean
}> = ({ children, nodeId, selected }) => {
  const theme = useTheme()
  return (
    <div
      style={{
        height: "100%",
        width: "100%",
        background: selected
          ? alpha(theme.palette.primary.light, 0.1)
          : undefined,
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
