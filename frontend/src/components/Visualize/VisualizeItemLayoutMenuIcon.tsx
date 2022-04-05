import React from 'react'
import { useDispatch, useSelector } from 'react-redux'
import IconButton from '@mui/material/IconButton'
import Menu from '@mui/material/Menu'
import MenuItem from '@mui/material/MenuItem'
import ListItemText from '@mui/material/ListItemText'
import ListItemIcon from '@mui/material/ListItemIcon'
import MoreVertIcon from '@mui/icons-material/MoreVert'
import DeleteIcon from '@mui/icons-material/Delete'
import AddIcon from '@mui/icons-material/Add'
import {
  selectVisualizeDataFilePath,
  selectVisualizeDataType,
} from 'store/slice/VisualizeItem/VisualizeItemSelectors'
import {
  deleteVisualizeItem,
  insertInitialItemToNextColumn,
} from 'store/slice/VisualizeItem/VisualizeItemSlice'
import { deleteDisplayItem } from 'store/slice/DisplayData/DisplayDataSlice'

export const VisualizeItemLayoutMenuIcon = React.memo<{ itemId: number }>(
  ({ itemId }) => {
    return <DisplayDataItemLayoutMenuIcon itemId={itemId} />
  },
)

const DisplayDataItemLayoutMenuIcon = React.memo<{
  itemId: number
}>(({ itemId }) => {
  const dispatch = useDispatch()
  const dataType = useSelector(selectVisualizeDataType(itemId))
  const filePath = useSelector(selectVisualizeDataFilePath(itemId))
  const onClickDeleteMenu = () => {
    dispatch(deleteVisualizeItem(itemId))
    // visualize Itemで同じpathのデータ個数を調べて、1だったら、displayも削除
    dispatch(deleteDisplayItem({ dataType, filePath }))
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
})

const PresentationalLayoutMenuIcon = React.memo<{
  onClickDeleteMenu: () => void
  onClickInsertMenu: () => void
}>(({ onClickDeleteMenu, onClickInsertMenu }) => {
  const [open, setOpen] = React.useState(false)
  const anchorRef = React.useRef<HTMLButtonElement>(null)
  const onClick: React.MouseEventHandler<HTMLButtonElement> = (e) => {
    e.stopPropagation() // 親divのonClickを反応させないため
    setOpen((prevOpen) => !prevOpen)
  }
  const onClose = () => {
    setOpen(false)
  }
  const onClickDeleteMenuFn: React.MouseEventHandler<HTMLLIElement> = (e) => {
    e.stopPropagation()
    onClickDeleteMenu()
    setOpen(false)
  }
  const onClickInsertMenuFn: React.MouseEventHandler<HTMLLIElement> = (e) => {
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
        <MenuItem onClick={onClickDeleteMenuFn}>
          <ListItemIcon>
            <DeleteIcon />
          </ListItemIcon>
          <ListItemText>Delete</ListItemText>
        </MenuItem>
        <MenuItem onClick={onClickInsertMenuFn}>
          <ListItemIcon>
            <AddIcon />
          </ListItemIcon>
          <ListItemText>Insert into next column</ListItemText>
        </MenuItem>
      </Menu>
    </>
  )
})
