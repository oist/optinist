import { FC } from "react"

import GitHubIcon from "@mui/icons-material/GitHub"
import MenuBookIcon from "@mui/icons-material/MenuBook"
import { Tooltip } from "@mui/material"
import IconButton from "@mui/material/IconButton"

const Tooltips: FC = () => {
  return (
    <>
      <Tooltip title="GitHub repository">
        <IconButton
          sx={{ mr: 1 }}
          href="https://github.com/oist/optinist"
          target="_blank"
        >
          <GitHubIcon />
        </IconButton>
      </Tooltip>
      <Tooltip title="Documentation">
        <IconButton
          href="https://optinist.readthedocs.io/en/latest/"
          target="_blank"
        >
          <MenuBookIcon />
        </IconButton>
      </Tooltip>
    </>
  )
}

export default Tooltips
