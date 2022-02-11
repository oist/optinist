import Dialog from '@mui/material/Dialog'
import { useDispatch, useSelector } from 'react-redux'
import TreeView from '@mui/lab/TreeView'
import TreeItem from '@mui/lab/TreeItem'
import FolderIcon from '@mui/icons-material/Folder'
import InsertDriveFileOutlinedIcon from '@mui/icons-material/InsertDriveFileOutlined'
import DialogActions from '@mui/material/DialogActions'
import DialogContent from '@mui/material/DialogContent'
import DialogTitle from '@mui/material/DialogTitle'
import Typography from '@mui/material/Typography'
import LinearProgress from '@mui/material/LinearProgress'
import Button from '@mui/material/Button'
import { useTheme } from '@mui/material/styles'
import React from 'react'

import {
  selectFilesIsLatest,
  selectFilesIsLoading,
  selectFilesTreeNodes,
} from 'store/slice/FilesTree/FilesTreeSelectors'
import { getFilesTree } from 'store/slice/FilesTree/FilesTreeAction'
import {
  FILE_TREE_TYPE,
  FILE_TREE_TYPE_SET,
  TreeNodeType,
} from 'store/slice/FilesTree/FilesTreeType'

type FileSelectDialogProps = {
  selectedFilePath: string
  onClickOk: (path: string) => void
  fileType?: FILE_TREE_TYPE
  title?: string
  open: boolean
  onClickCancel: () => void
  onClose?: () => void
}

export const FileSelectDialog = React.memo<FileSelectDialogProps>(
  ({
    open,
    selectedFilePath,
    onClickCancel,
    onClickOk,
    onClose,
    title,
    fileType = FILE_TREE_TYPE_SET.ALL,
  }) => {
    const [clickedFilePath, setClickedFilePath] =
      React.useState(selectedFilePath)
    const theme = useTheme()
    return (
      <Dialog open={open} onClose={onClose} fullWidth>
        <DialogTitle>{title ?? 'Select File'}</DialogTitle>
        <DialogContent dividers>
          <div
            style={{
              height: 300,
              overflow: 'auto',
              marginBottom: theme.spacing(1),
              border: '1px solid',
              padding: theme.spacing(1),
              borderColor: theme.palette.divider,
            }}
          >
            <FileTreeView
              onClickFile={setClickedFilePath}
              fileType={fileType}
            />
          </div>
          <Typography variant="subtitle1">Select File</Typography>
          <Typography variant="subtitle2">
            {!!clickedFilePath ? clickedFilePath : '---'}
          </Typography>
        </DialogContent>
        <DialogActions>
          <Button onClick={onClickCancel} variant="outlined">
            cancel
          </Button>
          <Button
            onClick={() => onClickOk(clickedFilePath)}
            color="primary"
            variant="outlined"
            autoFocus
          >
            OK
          </Button>
        </DialogActions>
      </Dialog>
    )
  },
)

const FileTreeView = React.memo<{
  onClickFile: (path: string) => void
  fileType: FILE_TREE_TYPE
}>(({ onClickFile, fileType }) => {
  const [tree, isLoading] = useFileTree(fileType)
  return (
    <div>
      {isLoading && <LinearProgress />}
      <TreeView>
        {tree?.map((node) => (
          <TreeNode node={node} onClickFile={onClickFile} />
        ))}
      </TreeView>
    </div>
  )
})

const TreeNode = React.memo<{
  node: TreeNodeType
  onClickFile: (path: string) => void
}>(({ node, onClickFile }) => {
  if (node.isDir) {
    return (
      <TreeItem
        icon={<FolderIcon htmlColor="skyblue" />}
        nodeId={node.path}
        label={node.name}
      >
        {node.nodes.map((childNode, i) => (
          <TreeNode node={childNode} key={i} onClickFile={onClickFile} />
        ))}
      </TreeItem>
    )
  } else {
    return (
      <TreeItem
        icon={<InsertDriveFileOutlinedIcon fontSize="small" />}
        nodeId={node.path}
        label={node.name}
        onClick={() => onClickFile(node.path)}
      />
    )
  }
})

function useFileTree(
  fileType: FILE_TREE_TYPE,
): [TreeNodeType[] | undefined, boolean] {
  const dispatch = useDispatch()
  const tree = useSelector(selectFilesTreeNodes(fileType))
  const isLatest = useSelector(selectFilesIsLatest(fileType))
  const isLoading = useSelector(selectFilesIsLoading(fileType))
  React.useEffect(() => {
    if (!isLatest && !isLoading) {
      dispatch(getFilesTree(fileType))
    }
  }, [isLatest, isLoading, fileType, dispatch])
  return [tree, isLoading]
}
