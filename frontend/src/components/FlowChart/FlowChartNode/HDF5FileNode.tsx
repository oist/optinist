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
} from '@material-ui/core'
import CloseOutlinedIcon from '@material-ui/icons/CloseOutlined'

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
  selectHDF5IsLatest,
  selectHDF5IsLoading,
  selectHDF5Nodes,
} from 'store/slice/HDF5/HDF5Selectors'
import { getHDF5Tree } from 'store/slice/HDF5/HDF5Action'
import { HDF5Tree, HDF5TreeDTO } from 'store/slice/HDF5/HDF5Type'
import { TreeItem, TreeView } from '@material-ui/lab'
import FolderIcon from '@material-ui/icons/Folder'
import InsertDriveFileOutlinedIcon from '@material-ui/icons/InsertDriveFileOutlined'

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
        id={toHandleId(nodeId, 'hdf5', 'HDF5Data')}
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
  const [tree, isLoading] = useHDF5Tree()
  console.log(tree)
  return (
    <div>
      {isLoading && <LinearProgress />}
      <TreeView>
        {tree?.map((node) => (
          <TreeNode node={node} />
        ))}
      </TreeView>
    </div>
  )
})

const TreeNode = React.memo<{
  node: HDF5TreeDTO
}>(({ node }) => {
  console.log(node)
  if (node.isDir) {
    // Directory
    return (
      <TreeItem
        icon={<FolderIcon htmlColor="skyblue" />}
        nodeId={node.name}
        label={node.name}
      >
        {node.nodes.map((childNode, i) => (
          <TreeNode node={childNode} key={i} />
        ))}
      </TreeItem>
    )
  } else {
    // File
    return (
      <TreeItem
        icon={<InsertDriveFileOutlinedIcon fontSize="small" />}
        nodeId={node.name}
        label={node.name + `, (shape=${node.shape})`}
      />
    )
  }
})

function useHDF5Tree(): [HDF5TreeDTO[] | undefined, boolean] {
  const dispatch = useDispatch()
  const tree = useSelector(selectHDF5Nodes())
  const isLatest = useSelector(selectHDF5IsLatest())
  const isLoading = useSelector(selectHDF5IsLoading())
  React.useEffect(() => {
    if (!isLatest && !isLoading) {
      dispatch(getHDF5Tree())
    }
  }, [isLatest, isLoading, dispatch])
  return [tree, isLoading]
}
