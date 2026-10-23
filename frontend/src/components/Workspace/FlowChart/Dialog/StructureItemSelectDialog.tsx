import React, { memo, useEffect, useState } from "react"
import { useDispatch, useSelector } from "react-redux"

import FolderIcon from "@mui/icons-material/Folder"
import InsertDriveFileOutlinedIcon from "@mui/icons-material/InsertDriveFileOutlined"
import { Box, Checkbox, Typography, Tooltip, Divider } from "@mui/material"
import Button from "@mui/material/Button"
import { CheckboxProps } from "@mui/material/Checkbox"
import Dialog from "@mui/material/Dialog"
import DialogActions from "@mui/material/DialogActions"
import DialogContent from "@mui/material/DialogContent"
import DialogTitle from "@mui/material/DialogTitle"
import LinearProgress from "@mui/material/LinearProgress"
import { useTheme } from "@mui/material/styles"
import { TreeItem } from "@mui/x-tree-view/TreeItem"
import { TreeView } from "@mui/x-tree-view/TreeView"

import {
  TreeNodeType,
  TreeDirType,
  TreeFileType,
  FileNodeConfig,
} from "components/Workspace/FlowChart/FlowChartNode/BaseStructuredFileNode"
import { NodeIdProps } from "store/slice/FlowElement/FlowElementType"
import { selectCurrentWorkspaceId } from "store/slice/Workspace/WorkspaceSelector"
import { AppDispatch } from "store/store"

type StructureItemSelectProps = {
  open: boolean
  setOpen: (value: boolean) => void
  config: FileNodeConfig
} & NodeIdProps

export const StructureItemSelectDialog = memo(
  function StructureItemSelectDialog({
    nodeId,
    open,
    setOpen,
    config,
  }: StructureItemSelectProps) {
    const dispatch = useDispatch<AppDispatch>()
    const [fileSelect, setFileSelect] = useState("")

    const structureFileName = useSelector(config.selectStructurePath(nodeId))

    const onClickOk = () => {
      dispatch(config.setStructurePath({ nodeId, path: fileSelect }))
      setOpen?.(false)
    }

    const onClickCancel = () => {
      setFileSelect("")
      setOpen?.(false)
    }

    return (
      <>
        <Typography className="selectFilePath" variant="caption">
          {structureFileName
            ? `↳ ${structureFileName}`
            : "No structure is selected."}
        </Typography>
        <Dialog open={open} onClose={() => setOpen(false)} fullWidth>
          <DialogTitle>Select Structure</DialogTitle>
          <StructureView
            nodeId={nodeId}
            fileSelect={fileSelect}
            setFileSelect={setFileSelect}
            config={config}
          />
          <DialogActions>
            <Button onClick={onClickCancel} color="primary" variant="outlined">
              cancel
            </Button>
            <Button onClick={onClickOk} variant="contained" autoFocus>
              OK
            </Button>
          </DialogActions>
        </Dialog>
      </>
    )
  },
)

const StructureView = memo(function StructureView({
  nodeId,
  fileSelect,
  setFileSelect,
  config,
}: NodeIdProps & {
  config: FileNodeConfig
  fileSelect: string
  setFileSelect: (value: string) => void
}) {
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
        <StructureTreeView
          nodeId={nodeId}
          fileSelect={fileSelect}
          setFileSelect={setFileSelect}
          config={config}
        />
      </div>
      <Typography>Selected Path</Typography>
      <Typography variant="subtitle2">{fileSelect || "---"}</Typography>
    </DialogContent>
  )
})

const StructureTreeView = memo(function StructureTreeView({
  nodeId,
  fileSelect,
  setFileSelect,
  config,
}: NodeIdProps & {
  config: FileNodeConfig
  fileSelect: string
  setFileSelect: (value: string) => void
}) {
  const [tree, isLoading] = useStructuredTree(nodeId, config)
  const [expanded, setExpanded] = useState<string[]>([])

  // Calculate paths to expand when fileSelect changes
  useEffect(() => {
    if (!fileSelect || !tree) return

    const pathsToExpand: string[] = []
    const findPathsToExpand = (
      nodes: TreeNodeType[],
      currentPath: string[] = [],
    ) => {
      for (const node of nodes) {
        if (node.isDir) {
          const dirNode = node as TreeDirType
          // Check if any child contains the selected file
          const containsSelectedFile = (n: TreeNodeType): boolean => {
            if (!n.isDir && n.path === fileSelect) return true
            if (n.isDir) {
              return (n as TreeDirType).nodes.some(containsSelectedFile)
            }
            return false
          }

          if (containsSelectedFile(node)) {
            pathsToExpand.push(node.path)
            findPathsToExpand(dirNode.nodes, [...currentPath, node.path])
          }
        }
      }
    }

    findPathsToExpand(tree)
    if (pathsToExpand.length > 0) {
      setExpanded(pathsToExpand)
    }
  }, [fileSelect, tree])

  const handleNodeToggle = (
    _event: React.SyntheticEvent,
    nodeIds: string[],
  ) => {
    setExpanded(nodeIds)
  }

  return (
    <div>
      {isLoading && <LinearProgress />}
      <Box display={"flex"} paddingBottom={1}>
        <Box flexGrow={4}>Structure</Box>
        <Box flexGrow={2}>Type</Box>
        <Box flexGrow={3}>Shape</Box>
        <Box flexGrow={2}>Nbytes</Box>
        <Box flexGrow={1}></Box>
      </Box>
      <Divider />
      <TreeView expanded={expanded} onNodeToggle={handleNodeToggle}>
        {tree?.map((node, i) => (
          <StructureTreeNode
            fileSelect={fileSelect}
            setFileSelect={setFileSelect}
            key={`${config.treeKeyPrefix}-${nodeId}-${i}`}
            node={node}
            nodeId={nodeId}
            config={config}
          />
        ))}
      </TreeView>
    </div>
  )
})

interface StructureTreeItemLabelProps {
  isFile: boolean
  shape: number[]
  type: string | null
  label: string
  nbytes?: string
  checkboxProps: CheckboxProps
}

const formatBytes = (nbytesStr?: string): string => {
  if (!nbytesStr) return ""

  const nbytes = parseInt(nbytesStr, 10)
  if (isNaN(nbytes)) return nbytesStr

  const BYTES_PER_KB = 1000
  const BYTES_PER_MB = BYTES_PER_KB * 1000

  if (nbytes < BYTES_PER_KB) {
    return `${nbytes} B`
  } else if (nbytes < BYTES_PER_MB) {
    const kb = nbytes / BYTES_PER_KB
    return `${kb.toFixed(1)} KB`
  } else {
    const mb = nbytes / BYTES_PER_MB
    return `${mb.toFixed(1)} MB`
  }
}

export const StructureTreeItemLabel = memo(function StructureTreeItemLabel({
  isFile = false,
  label,
  shape,
  type,
  nbytes,
  checkboxProps,
}: StructureTreeItemLabelProps) {
  return (
    <Box display="flex" alignItems="center" gap={2}>
      <Tooltip
        title={<span style={{ fontSize: 14 }}>{label}</span>}
        placement={"left"}
      >
        <Box
          width={isFile ? "25%" : "22%"}
          overflow={"hidden"}
          textOverflow={"ellipsis"}
        >
          {label}
        </Box>
      </Tooltip>
      <Box width={"15%"}>{type}</Box>
      <Box width={"25%"}>{shape ? `(${shape.join(", ")})` : ""}</Box>
      <Box width={"15%"}>{formatBytes(nbytes)}</Box>
      <Box>
        <Checkbox
          {...checkboxProps}
          disableRipple
          size="small"
          sx={{
            marginRight: "4px",
            padding: "2px",
          }}
        />
      </Box>
    </Box>
  )
})

interface StructureTreeNodeProps extends NodeIdProps {
  setFileSelect?: (value: string) => void
  fileSelect?: string
  node: TreeNodeType
  config: FileNodeConfig
}

const StructureTreeNode = memo(function TreeNode({
  node,
  nodeId,
  setFileSelect,
  fileSelect,
  config,
}: StructureTreeNodeProps) {
  const dispatch = useDispatch()
  const structureFileName = useSelector(config.selectStructurePath(nodeId))
  useEffect(() => {
    if (!structureFileName) return
    setFileSelect?.(structureFileName)
    //eslint-disable-next-line
  }, [dispatch, structureFileName])
  const onClickFile = (path: string) => {
    setFileSelect?.(path === fileSelect ? "" : path)
  }

  if (node.isDir) {
    // Directory
    return (
      <TreeItem
        icon={<FolderIcon htmlColor="skyblue" />}
        nodeId={node.path}
        label={node.name}
      >
        {(node as TreeDirType).nodes.map((childNode, i) => (
          <TreeNode
            setFileSelect={setFileSelect}
            fileSelect={fileSelect}
            node={childNode}
            key={i}
            nodeId={nodeId}
            config={config}
          />
        ))}
      </TreeItem>
    )
  } else {
    // File
    return (
      <TreeItem
        icon={<InsertDriveFileOutlinedIcon fontSize="small" />}
        nodeId={node.path}
        label={
          <StructureTreeItemLabel
            isFile={true}
            label={node.name}
            type={node.dataType || null}
            shape={(node as TreeFileType).shape || []}
            nbytes={(node as TreeFileType).nbytes}
            checkboxProps={{
              checked: fileSelect === node.path,
            }}
          />
        }
        onClick={() => onClickFile(node.path)}
      />
    )
  }
})

function useStructuredTree(
  nodeId: string,
  config: FileNodeConfig,
): [TreeNodeType[] | undefined, boolean] {
  const dispatch = useDispatch<AppDispatch>()
  const tree = useSelector(config.selectTree())
  const isLoading = useSelector(config.selectIsLoading())
  const filePathRaw = useSelector(config.selectFilePath(nodeId))
  const filePath = Array.isArray(filePathRaw) ? filePathRaw[0] : filePathRaw
  const workspaceId = useSelector(selectCurrentWorkspaceId)
  useEffect(() => {
    if (workspaceId && !isLoading && filePath) {
      dispatch(
        config.getTree({
          path: filePath as string,
          workspaceId: Number(workspaceId),
        }),
      )
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [workspaceId, filePath])
  return [tree, isLoading]
}
