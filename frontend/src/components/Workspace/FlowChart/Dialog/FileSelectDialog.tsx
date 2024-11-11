import {
  memo,
  SyntheticEvent,
  useContext,
  useEffect,
  useState,
  MouseEvent,
  useCallback,
} from "react"
import { useDispatch, useSelector } from "react-redux"

import AutorenewIcon from "@mui/icons-material/Autorenew"
import CloseIcon from "@mui/icons-material/Close"
import DeleteIcon from "@mui/icons-material/Delete"
import ErrorOutlineIcon from "@mui/icons-material/ErrorOutline"
import FolderIcon from "@mui/icons-material/Folder"
import InsertDriveFileOutlinedIcon from "@mui/icons-material/InsertDriveFileOutlined"
import { Divider, IconButton, Tooltip } from "@mui/material"
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
import { TreeItem } from "@mui/x-tree-view/TreeItem"
import { TreeView } from "@mui/x-tree-view/TreeView"

import { FILE_TREE_TYPE, FILE_TREE_TYPE_SET } from "api/files/Files"
import { DialogContext } from "components/Workspace/FlowChart/Dialog/DialogContext"
import {
  deleteFileTree,
  getFilesTree,
} from "store/slice/FilesTree/FilesTreeAction"
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
import { updateShape } from "store/slice/FileUploader/FileUploaderActions"
import { selectPipelineLatestUid } from "store/slice/Pipeline/PipelineSelectors"
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
  const { onOpenClearWorkflowIdDialog } = useContext(DialogContext)
  const currentWorkflowId = useSelector(selectPipelineLatestUid)
  const [selectedFilePath, setSelectedFilePath] = useState(initialFilePath)

  const onCancel = () => {
    setSelectedFilePath(initialFilePath) // 選択内容を反映させない
    onClickCancel()
  }
  const onOk = () => {
    if (currentWorkflowId != null) {
      onOpenClearWorkflowIdDialog({
        open: true,
        handleOk: () => {
          onClickOk(selectedFilePath)
        },
        handleCancel: () => onCancel(),
      })
    } else {
      onClickOk(selectedFilePath)
    }
  }
  const theme = useTheme()

  return (
    <Dialog open={open} onClose={onCancel} fullWidth>
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
        <Button onClick={onCancel} variant="outlined">
          cancel
        </Button>
        <Button onClick={onOk} variant="contained">
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
      {fileType === FILE_TREE_TYPE_SET.IMAGE ? (
        <>
          <Box sx={{ display: "flex" }}>
            <Typography sx={{ width: "50%" }}>Files</Typography>
            <Typography sx={{ minWidth: "35%" }} marginLeft={2}>
              Shapes
            </Typography>
            <Box></Box>
          </Box>
          <Divider />
        </>
      ) : null}
      <TreeView
        disableSelection={multiSelect}
        multiSelect={multiSelect}
        onNodeSelect={onNodeSelectHandler}
      >
        {tree?.map((node) => (
          <TreeNode
            fileType={fileType}
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
  fileType: FILE_TREE_TYPE
  node: TreeNodeType
  selectedFilePath: string[] | string
  multiSelect: boolean
  onCheckDir: (path: string, checked: boolean) => void
  onCheckFile: (path: string) => void
}

const TreeNode = memo(function TreeNode({
  fileType,
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
              isDir={node.isDir}
              fileType={fileType}
              shape={node.shape}
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
            fileType={fileType}
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
              isDir={node.isDir}
              fileType={fileType}
              shape={node.shape}
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
  fileType: FILE_TREE_TYPE
  shape: number[]
  label: string
  checkboxProps: CheckboxProps
  isDir?: boolean
}

const TreeItemLabel = memo(function TreeItemLabel({
  fileType,
  shape,
  label,
  isDir,
  checkboxProps,
}: TreeItemLabelProps) {
  const dispatch = useDispatch<AppDispatch>()
  const workspaceId = useSelector(selectCurrentWorkspaceId)
  const onUpdate = useCallback(
    (event: MouseEvent, fileName: string) => {
      if (!workspaceId) return
      event.stopPropagation()
      dispatch(updateShape({ workspaceId, fileName }))
    },
    [dispatch, workspaceId],
  )
  const onDelete = useCallback(
    (event: MouseEvent, fileName: string) => {
      if (!workspaceId) return
      event.stopPropagation()
      dispatch(deleteFileTree({ workspaceId, fileName, fileType }))
    },
    [dispatch, fileType, workspaceId],
  )
  return (
    <Box height={24} display="flex" alignItems="center">
      <Tooltip
        title={<span style={{ fontSize: 14 }}>{label}</span>}
        placement={"left-start"}
      >
        <Box
          sx={{
            width: "48%",
            textOverflow: "ellipsis",
            overflowX: "hidden",
            whiteSpace: "pre",
          }}
        >
          {label}
        </Box>
      </Tooltip>
      {fileType === FILE_TREE_TYPE_SET.IMAGE ? (
        <>
          <Box minWidth={175} marginLeft={2} alignItems={"center"}>
            {!isDir ? (
              !shape ? (
                <Tooltip
                  title={
                    <span style={{ fontSize: 14 }}>
                      parsing image shape failed
                    </span>
                  }
                  placement={"right"}
                >
                  <ErrorOutlineIcon color={"error"} />
                </Tooltip>
              ) : (
                `(${shape.join(", ")})`
              )
            ) : null}
          </Box>
        </>
      ) : null}
      <Box>
        <Checkbox
          {...checkboxProps}
          disableRipple
          size="small"
          sx={{
            marginRight: "4px",
            padding: "2px",
            minWidth: 24,
          }}
        />
      </Box>
      {!isDir ? (
        <IconButton
          sx={{ minWidth: 24 }}
          onClick={(event) => onUpdate(event, label)}
        >
          <AutorenewIcon />
        </IconButton>
      ) : (
        <Box width={24} marginRight={2} />
      )}
      <IconButton
        sx={{ minWidth: 24 }}
        onClick={(event) => onDelete(event, label)}
      >
        <DeleteIcon />
      </IconButton>
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
      {path ? (
        Array.isArray(path) ? (
          <ul style={{ padding: 0, margin: 0, listStyleType: "none" }}>
            {path.map((text) => (
              <li
                key={text}
                style={{
                  marginBottom: "8px",
                  listStyleType: "disc",
                  marginLeft: "24px",
                }}
              >
                <span
                  style={{
                    display: "flex",
                    justifyContent: "space-between",
                    width: "100%",
                  }}
                >
                  <span>{text}</span>
                  <IconButton style={{ padding: "0" }}>
                    <CloseIcon />
                  </IconButton>
                </span>
              </li>
            ))}
          </ul>
        ) : (
          path
        )
      ) : (
        "---"
      )}
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
