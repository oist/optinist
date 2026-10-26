import { FC, ReactNode } from "react"

import MenuIcon from "@mui/icons-material/Menu"
import MenuOpenIcon from "@mui/icons-material/MenuOpen"
import { Box, IconButton, Tooltip } from "@mui/material"
import { grey } from "@mui/material/colors"

import { DRAWER_WIDTH, CONTENT_HEIGHT } from "const/Layout"

interface LeftSidebarContainerProps {
  children: ReactNode
  isOpen: boolean
  onToggle: () => void
}

export const LeftSidebarContainer: FC<LeftSidebarContainerProps> = ({
  children,
  isOpen,
  onToggle,
}) => {
  return (
    <Box
      sx={{
        width: isOpen ? DRAWER_WIDTH : 56,
        height: CONTENT_HEIGHT,
        display: "flex",
        flexDirection: "column",
        marginRight: isOpen ? 3 : 0,
        borderRight: isOpen ? 1 : 0,
        borderColor: grey[300],
        overflow: "hidden",
        paddingTop: 0,
        paddingLeft: 1,
        transition: "width 0.3s ease-in-out, margin-right 0.3s ease-in-out",
        position: "relative",
      }}
    >
      <Box
        sx={{
          position: "absolute",
          top: 8,
          left: 8,
          opacity: isOpen ? 0 : 1,
          visibility: isOpen ? "hidden" : "visible",
          transition: "opacity 0.2s ease-in-out, visibility 0.2s ease-in-out",
          pointerEvents: isOpen ? "none" : "auto",
        }}
      >
        <Tooltip title="Open sidebar" placement="right">
          <IconButton
            onClick={onToggle}
            sx={{
              backgroundColor: "transparent",
              borderRadius: "6px",
              width: 36,
              height: 36,
              color: grey[700],
              "&:hover": {
                backgroundColor: grey[100],
                color: grey[900],
              },
            }}
          >
            <MenuIcon />
          </IconButton>
        </Tooltip>
      </Box>
      <Box
        sx={{
          display: "flex",
          flexDirection: "column",
          height: "100%",
          transform: isOpen ? "translateX(0)" : "translateX(-100%)",
          opacity: isOpen ? 1 : 0,
          transition: "transform 0.3s ease-in-out, opacity 0.3s ease-in-out",
          pointerEvents: isOpen ? "auto" : "none",
        }}
      >
        <Box
          sx={{
            display: "flex",
            justifyContent: "flex-start",
            alignItems: "center",
            paddingY: 0,
            borderBottom: `1px solid ${grey[200]}`,
            minHeight: 26,
          }}
        >
          <Tooltip title="Close sidebar">
            <IconButton
              onClick={onToggle}
              sx={{
                backgroundColor: "transparent",
                borderRadius: "6px",
                width: 36,
                height: 36,
                color: grey[700],
                "&:hover": {
                  backgroundColor: grey[100],
                  color: grey[900],
                },
              }}
            >
              <MenuOpenIcon />
            </IconButton>
          </Tooltip>
        </Box>
        <Box
          sx={{
            flex: 1,
            overflow: "auto",
          }}
        >
          {children}
        </Box>
      </Box>
    </Box>
  )
}
