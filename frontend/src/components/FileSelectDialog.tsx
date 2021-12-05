import Dialog from '@material-ui/core/Dialog'
import { useDispatch, useSelector } from 'react-redux'
import TreeView from '@material-ui/lab/TreeView'
import TreeItem from '@material-ui/lab/TreeItem'
import FolderIcon from '@material-ui/icons/Folder'
// import InsertDriveFileIcon from '@material-ui/icons/InsertDriveFile'
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
  filesTreeSelector,
} from 'store/slice/Files/FilesSelector'
import { getFiles } from 'store/slice/Files/FilesAction'
import { TreeNodeType } from 'store/slice/Files/FilesType'

type FileSelectDialogProps = {
  selectedFilePath: string
  onClickOk: (path: string) => void
  fileType?: string
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
    fileType,
  }) => {
    const [clickedFilePath, setClickedFilePath] =
      React.useState(selectedFilePath)
    const theme = useTheme()
    return (
      <Dialog open={open} onClose={onClose} fullWidth>
        <DialogTitle>{title ?? 'ファイルを選択'}</DialogTitle>
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
          <Typography variant="subtitle1">選択ファイル</Typography>
          <Typography variant="subtitle2">
            {!!clickedFilePath ? clickedFilePath : '---'}
          </Typography>
        </DialogContent>
        <DialogActions>
          <Button onClick={onClickCancel} variant="outlined">
            キャンセル
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
  fileType?: string
}>(({ onClickFile, fileType }) => {
  const [tree, isLoading] = useFileTree(fileType)
  return (
    <div>
      {isLoading && <LinearProgress />}
      <TreeView>
        {tree.map((node) => (
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

function useFileTree(fileType?: string): [TreeNodeType[], boolean] {
  const dispatch = useDispatch()
  const tree = useSelector(filesTreeSelector)
  const isLatest = useSelector(filesIsLatestSelector)
  const isLoading = useSelector(filesIsLoadingSelector)
  React.useEffect(() => {
    if (!isLatest && !isLoading) {
      dispatch(getFiles(fileType))
    }
  }, [isLatest, isLoading, fileType, dispatch])
  return [tree, isLoading]
}
