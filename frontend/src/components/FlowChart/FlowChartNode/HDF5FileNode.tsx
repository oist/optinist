import React, { CSSProperties } from 'react'
import { useDispatch, useSelector } from 'react-redux'
import { Handle, Position, NodeProps } from 'react-flow-renderer'

import { alpha, useTheme } from '@mui/material/styles'
import Dialog from '@mui/material/Dialog'
import TreeView from '@mui/lab/TreeView'
import TreeItem from '@mui/lab/TreeItem'
import FolderIcon from '@mui/icons-material/Folder'
import InsertDriveFileOutlinedIcon from '@mui/icons-material/InsertDriveFileOutlined'
import DialogActions from '@mui/material/DialogActions'
import DialogContent from '@mui/material/DialogContent'
import DialogTitle from '@mui/material/DialogTitle'
import LinearProgress from '@mui/material/LinearProgress'
import Button from '@mui/material/Button'

import { FILE_TYPE_SET } from 'store/slice/InputNode/InputNodeType'
import {
  selectHDF5InputNodeSelectedFilePath,
  selectInputNodeDefined,
  selectInputNodeHDF5Path,
} from 'store/slice/InputNode/InputNodeSelectors'
import { setInputNodeHDF5Path } from 'store/slice/InputNode/InputNodeSlice'
import { setInputNodeFilePath } from 'store/slice/InputNode/InputNodeActions'
import { toHandleId } from './FlowChartUtils'
import { FileSelect } from './FileSelect'
import { deleteFlowElementsById } from 'store/slice/FlowElement/FlowElementSlice'
import {
  selectHDF5IsLoading,
  selectHDF5Nodes,
} from 'store/slice/HDF5/HDF5Selectors'
import { getHDF5Tree } from 'store/slice/HDF5/HDF5Action'
import { HDF5TreeDTO } from 'store/slice/HDF5/HDF5Type'
import { Typography } from '@mui/material'

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
  const filePath = useSelector(selectHDF5InputNodeSelectedFilePath(nodeId))
  const onChangeFilePath = (path: string) => {
    dispatch(setInputNodeFilePath({ nodeId, filePath: path }))
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
      <button
        className="flowbutton"
        onClick={onClickDeleteIcon}
        style={{ color: 'black', position: 'absolute', top: -10, right: 10 }}
      >
        ×
      </button>
      <FileSelect
        onChangeFilePath={(path) => {
          if (!Array.isArray(path)) {
            onChangeFilePath(path)
          }
        }}
        fileType={FILE_TYPE_SET.HDF5}
        filePath={filePath ?? ''}
      />
      {filePath !== undefined && <ItemSelect nodeId={nodeId} />}
      <Handle
        type="source"
        position={Position.Right}
        id={toHandleId(nodeId, 'hdf5', 'HDF5Data')}
        style={sourceHandleStyle}
      />
    </div>
  )
})

const ItemSelect = React.memo<{
  nodeId: string
}>(({ nodeId }) => {
  const [open, setOpen] = React.useState(false)

  const structureFileName = useSelector(selectInputNodeHDF5Path(nodeId))

  return (
    <>
      <Button variant="outlined" size="small" onClick={() => setOpen(true)}>
        {'Structure'}
      </Button>
      <Typography className="selectFilePath" variant="caption">
        {!!structureFileName ? structureFileName : 'No structure is selected.'}
      </Typography>

      <Dialog open={open} onClose={() => setOpen(false)} fullWidth>
        <DialogTitle>{'Select File'}</DialogTitle>
        <Structure nodeId={nodeId} />
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

const Structure = React.memo<{
  nodeId: string
}>(({ nodeId }) => {
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
        <FileTreeView nodeId={nodeId} />
      </div>
    </DialogContent>
  )
})

const FileTreeView = React.memo<{
  nodeId: string
}>(({ nodeId }) => {
  const [tree, isLoading] = useHDF5Tree(nodeId)
  return (
    <div>
      {isLoading && <LinearProgress />}
      <TreeView>
        {tree?.map((node) => (
          <TreeNode node={node} nodeId={nodeId} />
        ))}
      </TreeView>
    </div>
  )
})

const TreeNode = React.memo<{
  node: HDF5TreeDTO
  nodeId: string
}>(({ node, nodeId }) => {
  const dispatch = useDispatch()

  const onClickFile = (path: string) => {
    dispatch(setInputNodeHDF5Path({ nodeId, path }))
  }

  if (node.isDir) {
    // Directory
    return (
      <TreeItem
        icon={<FolderIcon htmlColor="skyblue" />}
        nodeId={node.path}
        label={node.name}
      >
        {node.nodes.map((childNode, i) => (
          <TreeNode node={childNode} key={i} nodeId={nodeId} />
        ))}
      </TreeItem>
    )
  } else {
    // File
    return (
      <TreeItem
        icon={<InsertDriveFileOutlinedIcon fontSize="small" />}
        nodeId={node.path}
        label={node.name + `   (shape=${node.shape}, nbytes=${node.nbytes})`}
        onClick={() => onClickFile(node.path)}
      />
    )
  }
})

function useHDF5Tree(nodeId: string): [HDF5TreeDTO[] | undefined, boolean] {
  const dispatch = useDispatch()
  const tree = useSelector(selectHDF5Nodes())
  const isLoading = useSelector(selectHDF5IsLoading())
  const filePath = useSelector(selectHDF5InputNodeSelectedFilePath(nodeId))
  React.useEffect(() => {
    if (!isLoading && filePath) {
      dispatch(getHDF5Tree({ path: filePath }))
    }
  }, [isLoading, filePath, dispatch])
  return [tree, isLoading]
}
