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
import { isDirNodeByPath } from 'store/slice/FilesTree/FilesTreeUtils'

type FileSelectDialogProps = {
  initialFilePath: string[] | string
  onClickOk: (path: string[] | string) => void
  fileType?: FILE_TREE_TYPE
  title?: string
  open: boolean
  onClickCancel: () => void
  onClose?: () => void
  multiSelect: boolean
}

export const FileSelectDialog = React.memo<FileSelectDialogProps>(
  function FileSelectDialog({
    open,
    initialFilePath,
    onClickCancel,
    onClickOk,
    onClose,
    title,
    fileType = FILE_TREE_TYPE_SET.ALL,
    multiSelect,
  }) {
    const [selectedFilePath, setSelectedFilePath] =
      React.useState(initialFilePath)
    const onNodeSelect = (path: string[] | string) => {
      setSelectedFilePath(path)
    }
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
              onNodeSelect={onNodeSelect}
              multiSelect={multiSelect}
              fileType={fileType}
            />
          </div>
          <Typography variant="subtitle1">Select File</Typography>
          <FilePathSelectedListView path={selectedFilePath} />
        </DialogContent>
        <DialogActions>
          <Button onClick={onClickCancel} variant="outlined" color="inherit">
            cancel
          </Button>
          <Button
            onClick={() => onClickOk(selectedFilePath)}
            color="primary"
            variant="outlined"
          >
            OK
          </Button>
        </DialogActions>
      </Dialog>
    )
  },
)

const FilePathSelectedListView = React.memo<{ path: string | string[] }>(
  ({ path }) => {
    return (
      <Typography variant="subtitle2">
        {!!path
          ? Array.isArray(path)
            ? path.map((text) => <li>{text}</li>)
            : path
          : '---'}
      </Typography>
    )
  },
)

const FileTreeView = React.memo<{
  onNodeSelect: (path: string[] | string) => void
  multiSelect: boolean
  fileType: FILE_TREE_TYPE
}>(({ onNodeSelect, fileType, multiSelect }) => {
  const [tree, isLoading] = useFileTree(fileType)
  const onNodeSelectHandler = (
    event: React.SyntheticEvent,
    nodeIds: Array<string> | string,
  ) => {
    if (multiSelect) {
      onNodeSelect(nodeIds)
    } else {
      if (tree != null) {
        // multiSelectがfalseの場合、ディレクトリは選択しない
        const path = nodeIds as string
        if (!isDirNodeByPath(path, tree)) {
          onNodeSelect(path)
        }
      }
    }
  }
  return (
    <div>
      {isLoading && <LinearProgress />}
      <TreeView multiSelect={multiSelect} onNodeSelect={onNodeSelectHandler}>
        {tree?.map((node) => (
          <TreeNode node={node} />
        ))}
      </TreeView>
    </div>
  )
})

const TreeNode = React.memo<{
  node: TreeNodeType
}>(({ node }) => {
  if (node.isDir) {
    return (
      <TreeItem
        icon={<FolderIcon htmlColor="skyblue" />}
        nodeId={node.path}
        label={node.name}
      >
        {node.nodes.map((childNode, i) => (
          <TreeNode node={childNode} key={i} />
        ))}
      </TreeItem>
    )
  } else {
    return (
      <TreeItem
        icon={<InsertDriveFileOutlinedIcon fontSize="small" />}
        nodeId={node.path}
        label={node.name}
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
