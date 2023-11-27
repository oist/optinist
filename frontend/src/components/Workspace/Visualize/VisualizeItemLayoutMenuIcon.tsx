import { memo, useState, useRef, MouseEventHandler } from "react"
import { useDispatch, useSelector } from "react-redux"

import AddIcon from "@mui/icons-material/Add"
import DeleteIcon from "@mui/icons-material/Delete"
import MoreVertIcon from "@mui/icons-material/MoreVert"
import IconButton from "@mui/material/IconButton"
import ListItemIcon from "@mui/material/ListItemIcon"
import ListItemText from "@mui/material/ListItemText"
import Menu from "@mui/material/Menu"
import MenuItem from "@mui/material/MenuItem"

import { cancelRoiApi } from "api/outputs/Outputs"
import { deleteDisplayItem } from "store/slice/VisualizeItem/VisualizeItemActions"
import {
  selectDisplayDataIsSingle,
  selectRoiItemFilePath,
  selectVisualizeDataFilePath,
  selectVisualizeDataType,
} from "store/slice/VisualizeItem/VisualizeItemSelectors"
import { insertInitialItemToNextColumn } from "store/slice/VisualizeItem/VisualizeItemSlice"
import { selectCurrentWorkspaceId } from "store/slice/Workspace/WorkspaceSelector"

interface ItemIdProps {
  itemId: number
}

export const DisplayDataItemLayoutMenuIcon = memo(
  function DisplayDataItemLayoutMenuIcon({ itemId }: ItemIdProps) {
    const dispatch = useDispatch()
    const dataType = useSelector(selectVisualizeDataType(itemId))
    const filePath = useSelector(selectVisualizeDataFilePath(itemId))
    const isSingleData = useSelector(selectDisplayDataIsSingle(itemId))
    const roiFilePath = useSelector(selectRoiItemFilePath(0))
    const workspaceId = useSelector(selectCurrentWorkspaceId)

    const onClickDeleteMenu = () => {
      // visualize Itemで同じpathのデータ個数を調べて、1だったら、displayも削除
      dispatch(
        deleteDisplayItem(
          isSingleData && filePath != null && dataType != null
            ? { itemId, deleteData: true, filePath, dataType }
            : { itemId, deleteData: false },
        ),
      )
      if (!roiFilePath || !workspaceId) return
      cancelRoiApi(roiFilePath, workspaceId)
    }
    const onClickInsertMenu = () => {
      dispatch(insertInitialItemToNextColumn(itemId))
    }
    return (
      <PresentationalLayoutMenuIcon
        onClickDeleteMenu={onClickDeleteMenu}
        onClickInsertMenu={onClickInsertMenu}
      />
    )
  },
)

interface PresentationalLayoutMenuIconProps {
  onClickDeleteMenu: () => void
  onClickInsertMenu: () => void
}

const PresentationalLayoutMenuIcon = memo(
  function PresentationalLayoutMenuIcon({
    onClickDeleteMenu,
    onClickInsertMenu,
  }: PresentationalLayoutMenuIconProps) {
    const [open, setOpen] = useState(false)
    const anchorRef = useRef<HTMLButtonElement>(null)
    const onClick: MouseEventHandler<HTMLButtonElement> = (e) => {
      e.stopPropagation() // 親divのonClickを反応させないため
      setOpen((prevOpen) => !prevOpen)
    }
    const onClose = () => {
      setOpen(false)
    }
    const onClickDeleteMenuFn: MouseEventHandler<HTMLLIElement> = (e) => {
      e.stopPropagation()
      onClickDeleteMenu()
      setOpen(false)
    }
    const onClickInsertMenuFn: MouseEventHandler<HTMLLIElement> = (e) => {
      e.stopPropagation()
      onClickInsertMenu()
      setOpen(false)
    }

    return (
      <>
        <IconButton ref={anchorRef} onClick={onClick}>
          <MoreVertIcon />
        </IconButton>
        <Menu anchorEl={anchorRef.current} open={open} onClose={onClose}>
          <MenuItem onClick={onClickInsertMenuFn}>
            <ListItemIcon>
              <AddIcon />
            </ListItemIcon>
            <ListItemText>Insert into next column</ListItemText>
          </MenuItem>
          <MenuItem onClick={onClickDeleteMenuFn}>
            <ListItemIcon>
              <DeleteIcon />
            </ListItemIcon>
            <ListItemText>Delete</ListItemText>
          </MenuItem>
        </Menu>
      </>
    )
  },
)
