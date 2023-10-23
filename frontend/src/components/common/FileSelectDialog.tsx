import { memo, SyntheticEvent, useEffect, useState } from "react"
import { useDispatch, useSelector } from "react-redux"

import FolderIcon from "@mui/icons-material/Folder"
import InsertDriveFileOutlinedIcon from "@mui/icons-material/InsertDriveFileOutlined"
import Box from "@mui/material/Box"
import Button from "@mui/material/Button"
import Checkbox, { CheckboxProps } from "@mui/material/Checkbox"
import Dialog from "@mui/material/Dialog"
import DialogActions from "@mui/material/DialogActions"
import DialogContent from "@mui/material/DialogContent"
import DialogTitle from "@mui/material/DialogTitle"
import LinearProgress from "@mui/material/LinearProgress"
import { useTheme } from "@mui/material/styles"
import Typography from "@mui/material/Typography"
import { TreeView, TreeItem } from "@mui/x-tree-view"

import { FILE_TREE_TYPE, FILE_TREE_TYPE_SET } from "api/files/Files"
import { getFilesTree } from "store/slice/FilesTree/FilesTreeAction"
import {
  selectFilesIsLatest,
  selectFilesIsLoading,
  selectFilesTreeNodes,
} from "store/slice/FilesTree/FilesTreeSelectors"
import { TreeNodeType } from "store/slice/FilesTree/FilesTreeType"
import {
  getNodeByPath,
  isDirNodeByPath,
} from "store/slice/FilesTree/FilesTreeUtils"
import { selectCurrentWorkspaceId } from "store/slice/Workspace/WorkspaceSelector"
import { AppDispatch } from "store/store"

type FileSelectDialogProps = {
  initialFilePath: string[] | string
  onClickOk: (path: string[] | string) => void
  fileType?: FILE_TREE_TYPE
  title?: string
  open: boolean
  onClickCancel: () => void
  multiSelect: boolean
}

export const FileSelectDialog = memo(function FileSelectDialog({
  open,
  initialFilePath,
  onClickCancel,
  onClickOk,
  title,
  fileType = FILE_TREE_TYPE_SET.ALL,
  multiSelect,
}: FileSelectDialogProps) {
  useEffect(() => {
    setSelectedFilePath(initialFilePath)
  }, [initialFilePath])
  const [selectedFilePath, setSelectedFilePath] = useState(initialFilePath)
  const onCancel = () => {
    setSelectedFilePath(initialFilePath) // 選択内容を反映させない
    onClickCancel()
  }
  const onOk = () => {
    onClickOk(selectedFilePath)
  }
  const theme = useTheme()
  return (
    <Dialog open={open} fullWidth>
      <DialogTitle>{title ?? "Select File"}</DialogTitle>
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
          <FileTreeView
            setSelectedFilePath={setSelectedFilePath}
            multiSelect={multiSelect}
            fileType={fileType}
            selectedFilePath={selectedFilePath}
          />
        </div>
        <Typography variant="subtitle1">Select File</Typography>
        <FilePathSelectedListView path={selectedFilePath} />
      </DialogContent>
      <DialogActions>
        <Button onClick={onCancel} variant="outlined" color="inherit">
          cancel
        </Button>
        <Button onClick={onOk} color="primary" variant="outlined">
          OK
        </Button>
      </DialogActions>
    </Dialog>
  )
})

interface FileTreeViewProps {
  setSelectedFilePath: (path: string[] | string) => void
  selectedFilePath: string[] | string
  multiSelect: boolean
  fileType: FILE_TREE_TYPE
}

const FileTreeView = memo(function FileTreeView({
  setSelectedFilePath,
  selectedFilePath,
  fileType,
  multiSelect,
}: FileTreeViewProps) {
  const [tree, isLoading] = useFileTree(fileType)
  const onNodeSelectHandler = (
    event: SyntheticEvent,
    nodeIds: Array<string> | string,
  ) => {
    if (!multiSelect && tree != null) {
      // multiSelectがfalseの場合、ディレクトリは選択しない
      const path = nodeIds as string
      if (!isDirNodeByPath(path, tree)) {
        setSelectedFilePath(path)
      }
    }
  }
  // multiSelectでチェックボックスを使用する時用のハンドラ
  const onCheckFile = (path: string) => {
    if (Array.isArray(selectedFilePath)) {
      if (selectedFilePath.includes(path)) {
        setSelectedFilePath(
          selectedFilePath.filter((selectedPath) => path !== selectedPath),
        )
      } else {
        setSelectedFilePath(selectedFilePath.concat(path))
      }
    }
  }
  const onCheckDir = (path: string, checked: boolean) => {
    if (tree != null && Array.isArray(selectedFilePath)) {
      const node = getNodeByPath(path, tree)
      if (node != null && node.isDir) {
        const childrenFilePathList = node.nodes
          .filter((node) => !node.isDir)
          .map((node) => node.path)
        if (checked) {
          setSelectedFilePath(
            // concat時の重複を削除
            Array.from(new Set(selectedFilePath.concat(childrenFilePathList))),
          )
        } else {
          setSelectedFilePath(
            selectedFilePath.filter(
              (selectedPath) => !childrenFilePathList.includes(selectedPath),
            ),
          )
        }
      }
    }
  }
  return (
    <div>
      {isLoading && <LinearProgress />}
      <TreeView
        disableSelection={multiSelect}
        multiSelect={multiSelect}
        onNodeSelect={onNodeSelectHandler}
      >
        {tree?.map((node) => (
          <TreeNode
            key={node.name}
            node={node}
            selectedFilePath={selectedFilePath}
            multiSelect={multiSelect}
            onCheckDir={onCheckDir}
            onCheckFile={onCheckFile}
          />
        ))}
      </TreeView>
    </div>
  )
})

interface TreeNodeProps {
  node: TreeNodeType
  selectedFilePath: string[] | string
  multiSelect: boolean
  onCheckDir: (path: string, checked: boolean) => void
  onCheckFile: (path: string) => void
}

const TreeNode = memo(function TreeNode({
  node,
  selectedFilePath,
  multiSelect,
  onCheckDir,
  onCheckFile,
}: TreeNodeProps) {
  if (node.isDir) {
    const allChecked =
      Array.isArray(selectedFilePath) &&
      node.nodes
        .filter((node) => !node.isDir)
        .map((node) => node.path)
        .every((filePath) => selectedFilePath.includes(filePath))
    const allNotChecked =
      Array.isArray(selectedFilePath) &&
      node.nodes
        .filter((node) => !node.isDir)
        .map((node) => node.path)
        .every((filePath) => !selectedFilePath.includes(filePath))
    const indeterminate = !(allChecked || allNotChecked)
    return (
      <TreeItem
        icon={<FolderIcon htmlColor="skyblue" />}
        nodeId={node.path}
        label={
          multiSelect && node.nodes.filter((node) => !node.isDir).length > 0 ? (
            <TreeItemLabel
              label={node.name}
              checkboxProps={{
                indeterminate,
                checked: allChecked,
                onClick: (e) => {
                  e.stopPropagation() // on/offのクリックにつられてTreeを開閉させないようにする
                },
                onChange: (e) => onCheckDir(node.path, e.target.checked),
              }}
            />
          ) : (
            node.name
          )
        }
      >
        {node.nodes.map((childNode, i) => (
          <TreeNode
            node={childNode}
            selectedFilePath={selectedFilePath}
            key={i}
            multiSelect={multiSelect}
            onCheckDir={onCheckDir}
            onCheckFile={onCheckFile}
          />
        ))}
      </TreeItem>
    )
  } else {
    return (
      <TreeItem
        icon={<InsertDriveFileOutlinedIcon fontSize="small" />}
        nodeId={node.path}
        label={
          multiSelect ? (
            <TreeItemLabel
              label={node.name}
              checkboxProps={{
                checked:
                  Array.isArray(selectedFilePath) &&
                  selectedFilePath.includes(node.path),
                onChange: () => onCheckFile(node.path),
              }}
            />
          ) : (
            node.name
          )
        }
        onClick={() => onCheckFile(node.path)}
      />
    )
  }
})

interface TreeItemLabelProps {
  label: string
  checkboxProps: CheckboxProps
}

const TreeItemLabel = memo(function TreeItemLabel({
  label,
  checkboxProps,
}: TreeItemLabelProps) {
  return (
    <Box display="flex" alignItems="center">
      <Box flexGrow={1}>{label}</Box>
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

interface FilePathProps {
  path: string | string[]
}

const FilePathSelectedListView = memo(function FilePathSelectedListView({
  path,
}: FilePathProps) {
  return (
    <Typography variant="subtitle2">
      {path
        ? Array.isArray(path)
          ? path.map((text) => <li key={text}>{text}</li>)
          : path
        : "---"}
    </Typography>
  )
})

function useFileTree(
  fileType: FILE_TREE_TYPE,
): [TreeNodeType[] | undefined, boolean] {
  const dispatch = useDispatch<AppDispatch>()
  const tree = useSelector(selectFilesTreeNodes(fileType))
  const isLatest = useSelector(selectFilesIsLatest(fileType))
  const isLoading = useSelector(selectFilesIsLoading(fileType))
  const workspaceId = useSelector(selectCurrentWorkspaceId)
  useEffect(() => {
    if (workspaceId && !isLatest && !isLoading) {
      dispatch(getFilesTree({ workspaceId, fileType }))
    }
  }, [workspaceId, isLatest, isLoading, fileType, dispatch])
  return [tree, isLoading]
}
