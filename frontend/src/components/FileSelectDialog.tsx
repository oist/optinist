import Dialog from '@material-ui/core/Dialog'
import { useDispatch, useSelector } from 'react-redux'
import TreeView from '@material-ui/lab/TreeView'
import TreeItem from '@material-ui/lab/TreeItem'
import FolderIcon from '@material-ui/icons/Folder'
import InsertDriveFileOutlinedIcon from '@material-ui/icons/InsertDriveFileOutlined'
import DialogActions from '@material-ui/core/DialogActions'
import DialogContent from '@material-ui/core/DialogContent'
import DialogTitle from '@material-ui/core/DialogTitle'
import Typography from '@material-ui/core/Typography'
import LinearProgress from '@material-ui/core/LinearProgress'
import Button from '@material-ui/core/Button'
import { useTheme } from '@material-ui/core/styles'
import React from 'react'

import {
  filesIsLatestSelector,
  filesIsLoadingSelector,
  filesTreeNodesSelector,
} from 'store/slice/FilesTree/FilesTreeSelector'
import { getFilesTree } from 'store/slice/FilesTree/FilesTreeAction'
import {
  FILE_TYPE,
  FILE_TYPE_SET,
  TreeNodeType,
} from 'store/slice/FilesTree/FilesTreeType'

type FileSelectDialogProps = {
  selectedFilePath: string
  onClickOk: (path: string) => void
  fileType?: FILE_TYPE
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
    fileType = FILE_TYPE_SET.ALL,
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
  fileType: FILE_TYPE
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
  fileType: FILE_TYPE,
): [TreeNodeType[] | undefined, boolean] {
  const dispatch = useDispatch()
  const tree = useSelector(filesTreeNodesSelector(fileType))
  const isLatest = useSelector(filesIsLatestSelector(fileType))
  const isLoading = useSelector(filesIsLoadingSelector(fileType))
  React.useEffect(() => {
    if (!isLatest && !isLoading) {
      dispatch(getFilesTree(fileType))
    }
  }, [isLatest, isLoading, fileType, dispatch])
  return [tree, isLoading]
}
