import { memo, useEffect, useState } from "react"
import { useDispatch, useSelector } from "react-redux"
import { Handle, Position, NodeProps } from "reactflow"

import FolderIcon from "@mui/icons-material/Folder"
import InsertDriveFileOutlinedIcon from "@mui/icons-material/InsertDriveFileOutlined"
import { Typography } from "@mui/material"
import Button from "@mui/material/Button"
import Dialog from "@mui/material/Dialog"
import DialogActions from "@mui/material/DialogActions"
import DialogContent from "@mui/material/DialogContent"
import DialogTitle from "@mui/material/DialogTitle"
import LinearProgress from "@mui/material/LinearProgress"
import { useTheme } from "@mui/material/styles"
import { TreeItem } from "@mui/x-tree-view/TreeItem"
import { TreeView } from "@mui/x-tree-view/TreeView"

import { FileSelect } from "components/Workspace/FlowChart/FlowChartNode/FileSelect"
import { toHandleId } from "components/Workspace/FlowChart/FlowChartNode/FlowChartUtils"
import { NodeContainer } from "components/Workspace/FlowChart/FlowChartNode/NodeContainer"
import { HANDLE_STYLE } from "const/flowchart"
import { deleteFlowNodeById } from "store/slice/FlowElement/FlowElementSlice"
import { NodeIdProps } from "store/slice/FlowElement/FlowElementType"
import { getHDF5Tree } from "store/slice/HDF5/HDF5Action"
import {
  selectHDF5IsLoading,
  selectHDF5Nodes,
} from "store/slice/HDF5/HDF5Selectors"
import { HDF5TreeNodeType } from "store/slice/HDF5/HDF5Type"
import { setInputNodeFilePath } from "store/slice/InputNode/InputNodeActions"
import {
  selectHDF5InputNodeSelectedFilePath,
  selectInputNodeDefined,
  selectInputNodeHDF5Path,
} from "store/slice/InputNode/InputNodeSelectors"
import { setInputNodeHDF5Path } from "store/slice/InputNode/InputNodeSlice"
import { FILE_TYPE_SET } from "store/slice/InputNode/InputNodeType"
import { selectCurrentWorkspaceId } from "store/slice/Workspace/WorkspaceSelector"
import { AppDispatch } from "store/store"

type ItemSelectProps = {
  open: boolean
  setOpen: (value: boolean) => void
} & NodeIdProps

export const HDF5FileNode = memo(function HDF5FileNode(element: NodeProps) {
  const defined = useSelector(selectInputNodeDefined(element.id))
  if (defined) {
    return <HDF5FileNodeImple {...element} />
  } else {
    return null
  }
})

const HDF5FileNodeImple = memo(function HDF5FileNodeImple({
  id: nodeId,
  selected,
}: NodeProps) {
  const dispatch = useDispatch()
  const filePath = useSelector(selectHDF5InputNodeSelectedFilePath(nodeId))

  const [open, setOpen] = useState(false)
  const onChangeFilePath = (path: string) => {
    dispatch(setInputNodeFilePath({ nodeId, filePath: path }))
  }

  const onClickDeleteIcon = () => {
    dispatch(deleteFlowNodeById(nodeId))
  }

  return (
    <NodeContainer nodeId={nodeId} selected={selected}>
      <button
        className="flowbutton"
        onClick={onClickDeleteIcon}
        style={{ color: "black", position: "absolute", top: -10, right: 10 }}
      >
        ×
      </button>
      <FileSelect
        nodeId={nodeId}
        onChangeFilePath={(path) => {
          if (!Array.isArray(path)) {
            onChangeFilePath(path)
          }
        }}
        setOpen={setOpen}
        fileType={FILE_TYPE_SET.HDF5}
        filePath={filePath ?? ""}
      />
      {filePath !== undefined && (
        <ItemSelect open={open} setOpen={setOpen} nodeId={nodeId} />
      )}
      <Handle
        type="source"
        position={Position.Right}
        id={toHandleId(nodeId, "hdf5", "HDF5Data")}
        style={{ ...HANDLE_STYLE }}
      />
    </NodeContainer>
  )
})

const ItemSelect = memo(function ItemSelect({
  nodeId,
  open,
  setOpen,
}: ItemSelectProps) {
  const structureFileName = useSelector(selectInputNodeHDF5Path(nodeId))

  return (
    <>
      <Typography className="selectFilePath" variant="caption">
        {structureFileName ? structureFileName : "No structure is selected."}
      </Typography>
      <Dialog open={open} onClose={() => setOpen(false)} fullWidth>
        <DialogTitle>{"Select File"}</DialogTitle>
        <Structure nodeId={nodeId} />
        <DialogActions>
          <Button onClick={() => setOpen(false)} variant="outlined">
            cancel
          </Button>
          <Button
            onClick={() => setOpen(false)}
            color="primary"
            variant="contained"
            autoFocus
          >
            OK
          </Button>
        </DialogActions>
      </Dialog>
    </>
  )
})

const Structure = memo(function Structure({ nodeId }: NodeIdProps) {
  const theme = useTheme()
  return (
    <DialogContent dividers>
      <div
        style={{
          height: 300,
          overflow: "auto",
          marginBottom: theme.spacing(1),
          border: "1px solid",
          padding: theme.spacing(1),
          borderColor: theme.palette.divider,
        }}
      >
        <FileTreeView nodeId={nodeId} />
      </div>
    </DialogContent>
  )
})

const FileTreeView = memo(function FileTreeView({ nodeId }: NodeIdProps) {
  const [tree, isLoading] = useHDF5Tree(nodeId)
  return (
    <div>
      {isLoading && <LinearProgress />}
      <TreeView>
        {tree?.map((node, i) => (
          <TreeNode
            key={`hdf5tree-${nodeId}-${i}`}
            node={node}
            nodeId={nodeId}
          />
        ))}
      </TreeView>
    </div>
  )
})

interface TreeNodeProps extends NodeIdProps {
  node: HDF5TreeNodeType
}

const TreeNode = memo(function TreeNode({ node, nodeId }: TreeNodeProps) {
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

function useHDF5Tree(
  nodeId: string,
): [HDF5TreeNodeType[] | undefined, boolean] {
  const dispatch = useDispatch<AppDispatch>()
  const tree = useSelector(selectHDF5Nodes())
  const isLoading = useSelector(selectHDF5IsLoading())
  const filePath = useSelector(selectHDF5InputNodeSelectedFilePath(nodeId))
  const workspaceId = useSelector(selectCurrentWorkspaceId)
  useEffect(() => {
    if (workspaceId && !isLoading && filePath) {
      dispatch(getHDF5Tree({ path: filePath, workspaceId }))
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [workspaceId, filePath])
  return [tree, isLoading]
}
