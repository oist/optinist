import { FC, MouseEvent, useState } from "react"
import { useDispatch, useSelector } from "react-redux"

import { useSnackbar } from "notistack"

import { Addchart, GitHub, MenuBook, OpenInNew } from "@mui/icons-material"
import {
  ListItemIcon,
  ListItemText,
  Menu,
  MenuItem,
  Tooltip,
} from "@mui/material"
import IconButton from "@mui/material/IconButton"

import { getExperiments } from "store/slice/Experiments/ExperimentsActions"
import { reset } from "store/slice/VisualizeItem/VisualizeItemSlice"
import { importSampleData } from "store/slice/Workflow/WorkflowActions"
import { selectCurrentWorkspaceId } from "store/slice/Workspace/WorkspaceSelector"
import { AppDispatch } from "store/store"

const Tooltips: FC = () => {
  const [anchorEl, setAnchorEl] = useState<null | HTMLElement>(null)
  const menuId = "documentation-menu"
  const open = Boolean(anchorEl)
  const handleClickMenuIcon = (event: MouseEvent<HTMLElement>) => {
    setAnchorEl(event.currentTarget)
  }
  const handleClose = () => {
    setAnchorEl(null)
  }
  const handleGoToDocClick = () => {
    window.open("https://optinist.readthedocs.io/en/latest/", "_blank")
  }

  const dispatch: AppDispatch = useDispatch()
  const { enqueueSnackbar } = useSnackbar()
  const workspaceId = useSelector(selectCurrentWorkspaceId)
  const handleImportSampleDataClick = () => {
    if (typeof workspaceId === "number") {
      dispatch(importSampleData({ workspaceId }))
        .unwrap()
        .then(() => {
          enqueueSnackbar("Sample data import success", { variant: "success" })
          dispatch(reset())
          dispatch(getExperiments())
        })
        .catch(() => {
          enqueueSnackbar("Sample data import error", { variant: "error" })
        })
    }
  }

  return (
    <>
      <Tooltip title="GitHub repository">
        <IconButton href="https://github.com/oist/optinist" target="_blank">
          <GitHub />
        </IconButton>
      </Tooltip>
      <Tooltip title="Documentation">
        <IconButton onClick={handleClickMenuIcon}>
          <MenuBook />
        </IconButton>
      </Tooltip>
      <Menu
        id={menuId}
        anchorEl={anchorEl}
        open={open}
        onClose={handleClose}
        MenuListProps={{
          "aria-labelledby": menuId,
          role: "listbox",
        }}
      >
        <MenuItem onClick={handleGoToDocClick}>
          <ListItemIcon>
            <OpenInNew />
          </ListItemIcon>
          <ListItemText>Go to documentation page</ListItemText>
        </MenuItem>
        <MenuItem onClick={handleImportSampleDataClick}>
          <ListItemIcon>
            <Addchart />
          </ListItemIcon>
          <ListItemText>Import sample data</ListItemText>
        </MenuItem>
      </Menu>
    </>
  )
}

export default Tooltips
