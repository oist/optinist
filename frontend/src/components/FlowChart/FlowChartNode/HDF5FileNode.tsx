import React, { CSSProperties } from 'react'
import { useDispatch, useSelector } from 'react-redux'
import { Handle, Position, NodeProps } from 'react-flow-renderer'
import { alpha, useTheme } from '@material-ui/core/styles'
import {
  Button,
  Dialog,
  DialogActions,
  DialogContent,
  DialogTitle,
  IconButton,
  LinearProgress,
  Typography,
} from '@material-ui/core'
import CloseOutlinedIcon from '@material-ui/icons/CloseOutlined'
import FolderIcon from '@material-ui/icons/Folder'
import { TreeItem, TreeView } from '@material-ui/lab'
import InsertDriveFileOutlinedIcon from '@material-ui/icons/InsertDriveFileOutlined'

import { FILE_TYPE_SET } from 'store/slice/InputNode/InputNodeType'
import {
  selectInputNodeDefined,
  selectInputNodeSelectedFilePath,
} from 'store/slice/InputNode/InputNodeSelectors'
import { setInputNodeFilePath } from 'store/slice/InputNode/InputNodeSlice'
import { toHandleId } from './FlowChartUtils'
import { FileSelect } from './FileSelect'
import {
  deleteFlowElementsById,
  edifFlowElementsLabelById,
} from 'store/slice/FlowElement/FlowElementSlice'
import {
  FILE_TREE_TYPE,
  TreeNodeType,
} from 'store/slice/FilesTree/FilesTreeType'
import // selectHDF5IsLatest,
// selectHDF5IsLoading,
// selectHDF5TreeNodes,
'store/slice/FilesTree/FilesTreeSelectors'

const sourceHandleStyle: CSSProperties = {
  width: 8,
  height: 15,
  top: 15,
  border: '1px solid',
  borderColor: '#555',
  borderRadius: 0,
}

export const HDF5FileNode = React.memo<NodeProps>((element) => {
  const defined = useSelector(selectInputNodeDefined(element.id))
  if (defined) {
    return <HDF5FileNodeImple {...element} />
  } else {
    return null
  }
})

const HDF5FileNodeImple = React.memo<NodeProps>(({ id: nodeId, selected }) => {
  const dispatch = useDispatch()
  const filePath = useSelector(selectInputNodeSelectedFilePath(nodeId))
  const onChangeFilePath = (path: string) => {
    dispatch(setInputNodeFilePath({ nodeId, filePath: path }))
    const fileName = path.split('/').reverse()[0]
    dispatch(
      edifFlowElementsLabelById({
        nodeId,
        fileName,
      }),
    )
  }
  const theme = useTheme()

  const onClickDeleteIcon = () => {
    dispatch(deleteFlowElementsById(nodeId))
  }

  return (
    <div
      style={{
        height: '100%',
        width: '230px',
        background: selected
          ? alpha(theme.palette.primary.light, 0.1)
          : undefined,
      }}
    >
      <IconButton
        aria-label="delete"
        style={{ color: 'black', position: 'absolute', top: -20, right: -5 }}
        onClick={onClickDeleteIcon}
      >
        <CloseOutlinedIcon />
      </IconButton>
      <FileSelect
        onChangeFilePath={onChangeFilePath}
        fileType={FILE_TYPE_SET.HDF5}
        filePath={filePath ? filePath.split('/').reverse()[0] : ''}
      />
      <ItemSelect />
      <Handle
        type="source"
        position={Position.Right}
        id={toHandleId(nodeId, 'HDF5', 'HDF5Data')}
        style={sourceHandleStyle}
      />
    </div>
  )
})

const ItemSelect = React.memo(() => {
  const [open, setOpen] = React.useState(false)

  const onSelectFile = (selectedPath: string) => {
    // onChangeFilePath(selectedPath)
  }

  return (
    <>
      <Button variant="outlined" onClick={() => setOpen(true)}>
        {'Structure'}
      </Button>
      <Dialog open={open} onClose={() => setOpen(false)} fullWidth>
        <DialogTitle>{'Select File'}</DialogTitle>
        <Structure />
        <DialogActions>
          <Button onClick={() => setOpen(false)} variant="outlined">
            cancel
          </Button>
          <Button
            onClick={() => setOpen(false)}
            color="primary"
            variant="outlined"
            autoFocus
          >
            OK
          </Button>
        </DialogActions>
      </Dialog>
    </>
  )
})

const Structure = React.memo(() => {
  const theme = useTheme()
  return (
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
        <FileTreeView />
      </div>
    </DialogContent>
  )
})

const FileTreeView = React.memo<{}>(() => {
  // const [tree, isLoading] = useHDF5Tree()
  return <div></div>
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

function useHDF5Tree(): [] {
  const dispatch = useDispatch()
  // const tree = useSelector(selectHDF5TreeNodes())
  // const isLatest = useSelector(selectHDF5IsLatest())
  // const isLoading = useSelector(selectHDF5IsLoading())
  // React.useEffect(() => {
  //   if (!isLatest && !isLoading) {
  //     dispatch(getHDF5Tree())
  //   }
  // }, [isLatest, isLoading, dispatch])
  // return [tree, isLoading]
  return []
}
