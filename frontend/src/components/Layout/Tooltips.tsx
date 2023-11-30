import { FC } from "react"

import GitHubIcon from "@mui/icons-material/GitHub"
import MenuBookIcon from "@mui/icons-material/MenuBook"
import { Tooltip } from "@mui/material"
import IconButton from "@mui/material/IconButton"

import { ImportSampleDataButton } from "components/common/SampleData"

const Tooltips: FC = () => {
  return (
    <>
      <ImportSampleDataButton />
      <Tooltip title="GitHub repository">
        <IconButton
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
