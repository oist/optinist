import React, {
  memo,
  useContext,
  useEffect,
  useState,
  MouseEvent,
  useCallback,
  useRef,
  useMemo,
  createContext,
} from "react"
import { useDispatch, useSelector } from "react-redux"

import {
  DndContext,
  closestCenter,
  KeyboardSensor,
  PointerSensor,
  useSensor,
  useSensors,
  DragEndEvent,
} from "@dnd-kit/core"
import {
  arrayMove,
  SortableContext,
  sortableKeyboardCoordinates,
  verticalListSortingStrategy,
  useSortable,
} from "@dnd-kit/sortable"
import { CSS } from "@dnd-kit/utilities"
import AutorenewIcon from "@mui/icons-material/Autorenew"
import ClearIcon from "@mui/icons-material/Clear"
import CloseIcon from "@mui/icons-material/Close"
import DeleteIcon from "@mui/icons-material/Delete"
import DragIndicatorIcon from "@mui/icons-material/DragIndicator"
import ErrorOutlineIcon from "@mui/icons-material/ErrorOutline"
import FolderIcon from "@mui/icons-material/Folder"
import SearchIcon from "@mui/icons-material/Search"
import {
  Divider,
  IconButton,
  Tooltip,
  Chip,
  InputAdornment,
  TextField,
} from "@mui/material"
import Box from "@mui/material/Box"
import Button from "@mui/material/Button"
import Checkbox, { CheckboxProps } from "@mui/material/Checkbox"
import Dialog from "@mui/material/Dialog"
import DialogActions from "@mui/material/DialogActions"
import DialogContent from "@mui/material/DialogContent"
import DialogTitle from "@mui/material/DialogTitle"
import LinearProgress from "@mui/material/LinearProgress"
import { useTheme, styled } from "@mui/material/styles"
import Typography from "@mui/material/Typography"
import { TreeItem } from "@mui/x-tree-view/TreeItem"
import { TreeView } from "@mui/x-tree-view/TreeView"

import { FILE_TREE_TYPE, FILE_TREE_TYPE_SET } from "api/files/Files"
import { ConfirmDialog } from "components/common/ConfirmDialog"
import { DialogContext } from "components/Workspace/FlowChart/Dialog/DialogContext"
import { deleteFile, getFilesTree } from "store/slice/FilesTree/FilesTreeAction"
import {
  selectFilesIsLatest,
  selectFilesIsLoading,
  selectFilesTreeNodes,
} from "store/slice/FilesTree/FilesTreeSelectors"
import { TreeNodeType } from "store/slice/FilesTree/FilesTreeType"
import { getNodeByPath } from "store/slice/FilesTree/FilesTreeUtils"
import { updateShape } from "store/slice/FileUploader/FileUploaderActions"
import { selectPipelineLatestUid } from "store/slice/Pipeline/PipelineSelectors"
import { selectCurrentWorkspaceId } from "store/slice/Workspace/WorkspaceSelector"
import { AppDispatch } from "store/store"

// Constants
const COLUMN_MIN_WIDTH = 20
const COLUMN_MAX_WIDTH = 80
const COLUMN_DEFAULT_WIDTH = 50

// Context for file tree actions
type FileTreeActionsContextType = {
  onOpenDeleteDialog: (filePath: string, fileName: string) => void
}
const FileTreeActionsContext = createContext<FileTreeActionsContextType | null>(
  null,
)

// Common style objects for reuse
const commonStyles = {
  ellipsis: {
    overflow: "hidden",
    textOverflow: "ellipsis",
    whiteSpace: "nowrap",
  },
  flexCenter: {
    display: "flex",
    alignItems: "center",
  },
  iconButton: {
    minWidth: 24,
    padding: "2px",
  },
} as const

// Styled components (5+ attributes or commonly used)
const StyledColumnResizer = styled(Box)({
  width: "10px",
  cursor: "col-resize",
  position: "absolute",
  top: 0,
  bottom: "6px",
  zIndex: 10,
  display: "flex",
  alignItems: "center",
  justifyContent: "center",
  "&::before": {
    content: "''",
    width: "1px",
    height: "100%",
    backgroundColor: "rgba(0, 0, 0, 0.2)",
    display: "block",
  },
  "&:hover::before": {
    backgroundColor: "rgba(0, 0, 0, 0.4)",
    width: "2px",
  },
})

const StyledCheckbox = styled(Checkbox)({
  marginRight: "4px",
  padding: "2px",
  minWidth: 24,
})

// Custom Hooks
const useFilteredTree = (
  tree: TreeNodeType[] | null | undefined,
  filterText: string,
) => {
  const matchesFilter = useCallback((path: string, filter: string): boolean => {
    if (!filter) return true
    const pattern = filter
      .replace(/[.+?^${}()|[\]\\]/g, "\\$&")
      .replace(/\*/g, ".*")
    const regex = new RegExp(pattern, "i")
    return regex.test(path)
  }, [])

  return useMemo(() => {
    if (!tree)
      return { filteredTree: null, totalFileCount: 0, filteredFileCount: 0 }

    let totalCount = 0
    let filteredCount = 0

    const filterNodes = (nodes: TreeNodeType[]): TreeNodeType[] => {
      return nodes
        .map((node) => {
          if (node.isDir) {
            const filteredChildren = filterNodes(node.nodes)
            if (filteredChildren.length > 0) {
              return { ...node, nodes: filteredChildren }
            }
            return null
          } else {
            totalCount++
            if (matchesFilter(node.path, filterText)) {
              filteredCount++
              return node
            }
            return null
          }
        })
        .filter((node): node is TreeNodeType => node !== null)
    }

    const filtered = filterNodes(tree)
    return {
      filteredTree: filtered,
      totalFileCount: totalCount,
      filteredFileCount: filteredCount,
    }
  }, [tree, filterText, matchesFilter])
}

const useColumnResize = (initialWidth = COLUMN_DEFAULT_WIDTH) => {
  const [columnWidth, setColumnWidth] = useState(initialWidth)
  const [isDragging, setIsDragging] = useState(false)
  const containerRef = useRef<HTMLDivElement>(null)

  const handleMouseDown = useCallback((e: React.MouseEvent) => {
    e.preventDefault()
    setIsDragging(true)
  }, [])

  const handleMouseUp = useCallback(() => {
    setIsDragging(false)
  }, [])

  const handleMouseMove = useCallback(
    (e: Event) => {
      if (!isDragging || !containerRef.current) return
      const mouseEvent = e as unknown as MouseEvent
      const rect = containerRef.current.getBoundingClientRect()
      const newWidth = ((mouseEvent.clientX - rect.left) / rect.width) * 100
      setColumnWidth(
        Math.max(COLUMN_MIN_WIDTH, Math.min(COLUMN_MAX_WIDTH, newWidth)),
      )
    },
    [isDragging],
  )

  useEffect(() => {
    if (isDragging) {
      document.addEventListener("mousemove", handleMouseMove)
      document.addEventListener("mouseup", handleMouseUp)
      document.body.style.cursor = "col-resize"
      document.body.style.userSelect = "none"
      return () => {
        document.removeEventListener("mousemove", handleMouseMove)
        document.removeEventListener("mouseup", handleMouseUp)
        document.body.style.cursor = ""
        document.body.style.userSelect = ""
      }
    }
    return undefined
  }, [isDragging, handleMouseMove, handleMouseUp])

  return { columnWidth, handleMouseDown, containerRef }
}

const useFileSelection = (
  selectedFilePath: string[] | string,
  setSelectedFilePath: (path: string[] | string) => void,
  filteredTree: TreeNodeType[] | null | undefined,
) => {
  const getAllFiles = useCallback((): string[] => {
    if (!filteredTree) return []
    const files: string[] = []
    const collectFiles = (nodes: TreeNodeType[]) => {
      nodes.forEach((node) => {
        if (!node.isDir) {
          files.push(node.path)
        } else if (node.nodes) {
          collectFiles(node.nodes)
        }
      })
    }
    collectFiles(filteredTree)
    return files
  }, [filteredTree])

  const checkStatus = useMemo(() => {
    if (!Array.isArray(selectedFilePath) || !filteredTree)
      return { allChecked: false, someChecked: false }

    const allFiles = getAllFiles()
    if (allFiles.length === 0) return { allChecked: false, someChecked: false }

    const checkedCount = allFiles.filter((file) =>
      selectedFilePath.includes(file),
    ).length

    return {
      allChecked: checkedCount === allFiles.length,
      someChecked: checkedCount > 0 && checkedCount < allFiles.length,
    }
  }, [selectedFilePath, filteredTree, getAllFiles])

  const handleCheckAll = useCallback(
    (checked: boolean) => {
      if (checked) {
        const allFiles = getAllFiles()
        setSelectedFilePath(allFiles)
      } else {
        setSelectedFilePath([])
      }
    },
    [getAllFiles, setSelectedFilePath],
  )

  return { checkStatus, handleCheckAll, getAllFiles }
}

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
        <Box display="flex" alignItems="center" gap={1}>
          <Typography variant="subtitle1">
            {Array.isArray(selectedFilePath) && selectedFilePath.length === 0
              ? "No Selected File"
              : "Selected Files"}
          </Typography>
          {Array.isArray(selectedFilePath) && selectedFilePath.length > 0 && (
            <Chip
              label={selectedFilePath.length}
              size="small"
              color="primary"
              variant="outlined"
              sx={{ fontSize: "0.75rem", height: "20px", fontWeight: "bold" }}
            />
          )}
        </Box>
        <FilePathSelectedListView
          path={selectedFilePath}
          setSelectedFilePath={setSelectedFilePath}
        />
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
  const [filterText, setFilterText] = useState("") // Filter text state
  const [deleteDialogOpen, setDeleteDialogOpen] = useState(false)
  const [deleteTargetFile, setDeleteTargetFile] = useState<{
    path: string
    name: string
  } | null>(null)

  // Use custom hooks
  const { columnWidth, handleMouseDown, containerRef } = useColumnResize()
  const { filteredTree, totalFileCount, filteredFileCount } = useFilteredTree(
    tree,
    filterText,
  )
  const { checkStatus, handleCheckAll } = useFileSelection(
    selectedFilePath,
    setSelectedFilePath,
    filteredTree,
  )

  // Helper function to check if a file exists in the tree
  const isFileInTree = useCallback(
    (path: string, tree: TreeNodeType[] | null): boolean => {
      if (!tree) return false
      const checkNode = (node: TreeNodeType): boolean => {
        if (node.path === path) return true
        if (node.isDir && node.nodes) {
          return node.nodes.some(checkNode)
        }
        return false
      }
      return tree.some(checkNode)
    },
    [],
  )

  // Effect to remove selected file if it's not in the tree
  useEffect(() => {
    if (!tree) return

    const filterValidPaths = (paths: string | string[]) => {
      if (Array.isArray(paths)) {
        return paths.filter((p) => isFileInTree(p, tree))
      }
      return isFileInTree(paths, tree) ? paths : ""
    }

    const filtered = filterValidPaths(selectedFilePath)
    if (JSON.stringify(filtered) !== JSON.stringify(selectedFilePath)) {
      setSelectedFilePath(filtered)
    }
  }, [tree, isFileInTree, selectedFilePath, setSelectedFilePath])

  // File selection handlers
  const handleFileToggle = useCallback(
    (path: string) => {
      if (!Array.isArray(selectedFilePath)) {
        setSelectedFilePath(path === selectedFilePath ? "" : path)
        return
      }
      if (Array.isArray(selectedFilePath)) {
        setSelectedFilePath(
          selectedFilePath.includes(path)
            ? selectedFilePath.filter((p: string) => p !== path)
            : [...selectedFilePath, path],
        )
      }
    },
    [selectedFilePath, setSelectedFilePath],
  )

  const onCheckDir = useCallback(
    (path: string, checked: boolean) => {
      if (tree != null && Array.isArray(selectedFilePath)) {
        const node = getNodeByPath(path, tree)
        if (node != null && node.isDir) {
          const childrenFilePathList = node.nodes
            .filter((node) => !node.isDir)
            .map((node) => node.path)
          if (checked) {
            setSelectedFilePath(
              // concat時の重複を削除
              Array.from(
                new Set(selectedFilePath.concat(childrenFilePathList)),
              ),
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
    },
    [tree, selectedFilePath, setSelectedFilePath],
  )

  // Delete functionality
  const dispatch = useDispatch<AppDispatch>()
  const workspaceId = useSelector(selectCurrentWorkspaceId)

  const onOpenDeleteDialog = useCallback(
    (filePath: string, fileName: string) => {
      setDeleteTargetFile({ path: filePath, name: fileName })
      setDeleteDialogOpen(true)
    },
    [],
  )

  const onDeleteFile = useCallback(() => {
    if (!workspaceId || !deleteTargetFile) return

    // Remove the file from selectedFile state
    if (Array.isArray(selectedFilePath)) {
      setSelectedFilePath(
        selectedFilePath.filter((file) => file !== deleteTargetFile.path),
      )
    }

    dispatch(
      deleteFile({ workspaceId, fileName: deleteTargetFile.path, fileType }),
    )
    setDeleteDialogOpen(false)
    setDeleteTargetFile(null)
  }, [
    dispatch,
    fileType,
    selectedFilePath,
    setSelectedFilePath,
    workspaceId,
    deleteTargetFile,
  ])

  // All column resize and filter logic is now handled by custom hooks

  return (
    <div
      ref={containerRef}
      style={{ display: "flex", flexDirection: "column", height: "100%" }}
    >
      {isLoading && <LinearProgress />}
      {/* Filter input */}
      <Box sx={{ px: 1, pt: 0.75, pb: 2 }}>
        <TextField
          fullWidth
          size="small"
          placeholder="Filter... (* as wildcard)"
          value={filterText}
          onChange={(e) => setFilterText(e.target.value)}
          InputProps={{
            startAdornment: (
              <InputAdornment position="start">
                <SearchIcon fontSize="small" />
              </InputAdornment>
            ),
            endAdornment: filterText && (
              <InputAdornment position="end">
                <IconButton
                  size="small"
                  onClick={() => setFilterText("")}
                  edge="end"
                >
                  <ClearIcon fontSize="small" />
                </IconButton>
              </InputAdornment>
            ),
          }}
        />
        {/* File count info */}
        {filterText && (
          <Box
            sx={{ mt: 0.5, display: "flex", alignItems: "center", gap: 0.5 }}
          >
            <Typography variant="caption" color="text.secondary">
              Showing
            </Typography>
            <Chip
              label={filteredFileCount}
              size="small"
              sx={{ height: "18px", fontSize: "0.7rem" }}
              color={filteredFileCount === 0 ? "default" : "primary"}
              variant="outlined"
            />
            <Typography variant="caption" color="text.secondary">
              of {totalFileCount} files
            </Typography>
          </Box>
        )}
      </Box>
      <>
        <Box
          sx={{
            display: "flex",
            position: "sticky",
            top: 0,
            backgroundColor: "background.paper",
            zIndex: 5,
            userSelect: "none",
            paddingBottom: "6px",
            borderTop: "1px solid",
            borderTopColor: "divider",
            px: 1,
            pt: "6px",
          }}
        >
          <Box
            sx={{
              width: `${columnWidth}%`,
              display: "flex",
              alignItems: "center",
            }}
          >
            {multiSelect && (
              <StyledCheckbox
                checked={checkStatus.allChecked}
                indeterminate={checkStatus.someChecked}
                onChange={(e) => handleCheckAll(e.target.checked)}
                size="small"
                disableRipple
              />
            )}
            <Typography sx={{ ...commonStyles.ellipsis }}>Files</Typography>
          </Box>
          <StyledColumnResizer
            sx={{ left: `calc(${columnWidth}% - 5px)` }}
            onMouseDown={handleMouseDown}
          />
          {fileType === FILE_TREE_TYPE_SET.IMAGE && (
            <Typography
              sx={{ flex: 1, ...commonStyles.ellipsis }}
              marginLeft={2}
            >
              Shapes
            </Typography>
          )}
        </Box>
        <Divider />
      </>
      <FileTreeActionsContext.Provider value={{ onOpenDeleteDialog }}>
        <Box sx={{ flex: 1, overflow: "auto", px: 1 }}>
          <TreeView disableSelection={multiSelect} multiSelect={multiSelect}>
            {filteredTree?.map((node) => (
              <FileTreeNode
                fileType={fileType}
                key={node.name}
                node={node}
                selectedFilePath={selectedFilePath}
                multiSelect={multiSelect}
                onCheckDir={onCheckDir}
                onCheckFile={handleFileToggle}
                setSelectedFilePath={setSelectedFilePath}
                columnWidth={columnWidth}
              />
            ))}
          </TreeView>
        </Box>
      </FileTreeActionsContext.Provider>
      <ConfirmDialog
        open={deleteDialogOpen}
        setOpen={setDeleteDialogOpen}
        onConfirm={onDeleteFile}
        title="Are you sure you want to delete this item?"
        content={deleteTargetFile?.name || ""}
        confirmLabel="delete"
        iconType="warning"
      />
    </div>
  )
})

interface FileTreeNodeProps {
  fileType: FILE_TREE_TYPE
  node: TreeNodeType
  selectedFilePath: string[] | string
  multiSelect: boolean
  onCheckDir: (path: string, checked: boolean) => void
  onCheckFile: (path: string) => void
  setSelectedFilePath: (path: string[] | string) => void
  columnWidth?: number
}

const FileTreeNode = memo(function FileTreeNode({
  fileType,
  node,
  selectedFilePath,
  multiSelect,
  onCheckDir,
  onCheckFile,
  setSelectedFilePath,
  columnWidth = COLUMN_DEFAULT_WIDTH,
}: FileTreeNodeProps) {
  // Common tree item props
  const treeItemProps = {
    nodeId: node.path,
    sx: {
      "& .MuiTreeItem-iconContainer": {
        width: 0,
        minWidth: 0,
      },
    },
  }

  // Calculate checkbox props for directories
  const getDirCheckboxProps = useCallback(() => {
    if (!node.isDir || !multiSelect) return undefined

    const fileNodes = node.nodes.filter((n) => !n.isDir)
    if (fileNodes.length === 0) return undefined

    const filePaths = fileNodes.map((n) => n.path)
    const allChecked =
      Array.isArray(selectedFilePath) &&
      filePaths.every((path) => selectedFilePath.includes(path))
    const someChecked =
      Array.isArray(selectedFilePath) &&
      filePaths.some((path) => selectedFilePath.includes(path))

    return {
      indeterminate: someChecked && !allChecked,
      checked: allChecked,
      onClick: (e: React.MouseEvent) => e.stopPropagation(),
      onChange: (e: React.ChangeEvent<HTMLInputElement>) =>
        onCheckDir(node.path, e.target.checked),
    }
  }, [node, multiSelect, selectedFilePath, onCheckDir])

  // Calculate checkbox props for files
  const getFileCheckboxProps = useCallback(() => {
    const isChecked = multiSelect
      ? Array.isArray(selectedFilePath) && selectedFilePath.includes(node.path)
      : selectedFilePath === node.path

    return {
      checked: isChecked,
      onChange: (e: React.ChangeEvent<HTMLInputElement>) => {
        e.stopPropagation()
        onCheckFile(node.path)
      },
    }
  }, [multiSelect, selectedFilePath, node.path, onCheckFile])

  // Common label props
  const labelProps = {
    multiSelect,
    isDir: node.isDir,
    fileType,
    shape: node.shape,
    label: node.name,
    columnWidth,
    filePath: node.path,
    icon: node.isDir ? <FolderIcon htmlColor="skyblue" /> : undefined,
    checkboxProps: node.isDir ? getDirCheckboxProps() : getFileCheckboxProps(),
  }

  return (
    <TreeItem {...treeItemProps} label={<FileTreeItemLabel {...labelProps} />}>
      {node.isDir &&
        node.nodes.map((childNode, i) => (
          <FileTreeNode
            key={i}
            fileType={fileType}
            node={childNode}
            selectedFilePath={selectedFilePath}
            multiSelect={multiSelect}
            onCheckDir={onCheckDir}
            onCheckFile={onCheckFile}
            setSelectedFilePath={setSelectedFilePath}
            columnWidth={columnWidth}
          />
        ))}
    </TreeItem>
  )
})

interface FileTreeItemLabelProps {
  fileType: FILE_TREE_TYPE
  shape: number[]
  label: string
  checkboxProps?: CheckboxProps
  isDir?: boolean
  icon?: React.ReactNode
  multiSelect: boolean
  columnWidth?: number
  filePath: string
}

export const FileTreeItemLabel = memo(function FileTreeItemLabel({
  fileType,
  shape,
  label,
  isDir,
  checkboxProps,
  icon,
  multiSelect,
  columnWidth = COLUMN_DEFAULT_WIDTH,
  filePath,
}: FileTreeItemLabelProps) {
  const dispatch = useDispatch<AppDispatch>()
  const workspaceId = useSelector(selectCurrentWorkspaceId)
  const fileTreeActions = useContext(FileTreeActionsContext)

  // Consolidated file actions
  const fileActions = useMemo(
    () => ({
      updateShape: (event: MouseEvent) => {
        if (!workspaceId) return
        event.stopPropagation()
        dispatch(updateShape({ workspaceId, fileName: filePath }))
      },
      deleteFile: (event: MouseEvent) => {
        if (!workspaceId || !fileTreeActions) return
        event.stopPropagation()
        fileTreeActions.onOpenDeleteDialog(filePath, label)
      },
    }),
    [dispatch, workspaceId, fileTreeActions, filePath, label],
  )

  return (
    <Box
      sx={{
        "&:hover": {
          backgroundColor: "transparent",
        },
        cursor: "default",
      }}
    >
      <Box
        height={24}
        display="flex"
        justifyContent="space-between"
        alignItems="center"
      >
        <Tooltip
          title={<span style={{ fontSize: 14 }}>{label}</span>}
          placement={"left-start"}
        >
          <Box
            sx={{
              width: `${columnWidth}%`,
              display: "flex",
              alignItems: "center",
            }}
          >
            {checkboxProps && (
              <StyledCheckbox {...checkboxProps} size="small" disableRipple />
            )}
            {icon && (
              <Box sx={{ display: "flex", alignItems: "center", mr: 0.5 }}>
                {icon}
              </Box>
            )}
            <Box sx={{ ...commonStyles.ellipsis, whiteSpace: "pre", flex: 1 }}>
              {label}
            </Box>
          </Box>
        </Tooltip>
        {fileType === FILE_TREE_TYPE_SET.IMAGE ? (
          <>
            <Box
              flex={1}
              marginLeft={2}
              alignItems="center"
              sx={commonStyles.ellipsis}
            >
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
                  <span>{`(${shape.join(", ")})`}</span>
                )
              ) : null}
            </Box>
          </>
        ) : null}
        <Box display="flex" alignItems="center">
          {!isDir && multiSelect && (
            <IconButton
              sx={commonStyles.iconButton}
              onClick={fileActions.updateShape}
            >
              <AutorenewIcon />
            </IconButton>
          )}
          <IconButton
            sx={commonStyles.iconButton}
            color="error"
            onClick={fileActions.deleteFile}
            disabled={checkboxProps?.checked}
            data-testid="DeleteIconBtn"
          >
            <DeleteIcon />
          </IconButton>
        </Box>
      </Box>
    </Box>
  )
})

interface SortableItemProps {
  id: string
  text: string
  onRemove: (text: string) => void
}

const SortableItem = memo(function SortableItem({
  id,
  text,
  onRemove,
}: SortableItemProps) {
  const {
    attributes,
    listeners,
    setNodeRef,
    transform,
    transition,
    isDragging,
  } = useSortable({ id })

  const style = {
    transform: CSS.Transform.toString(transform),
    transition,
    opacity: isDragging ? 0.5 : 1,
  }

  return (
    <li
      ref={setNodeRef}
      style={{
        ...style,
        marginBottom: "3px",
        listStyleType: "none",
        marginLeft: "8px",
        cursor: isDragging ? "grabbing" : "grab",
      }}
    >
      <span
        style={{
          display: "flex",
          justifyContent: "space-between",
          alignItems: "center",
          width: "100%",
        }}
      >
        <span style={{ display: "flex", alignItems: "center", flex: 1 }}>
          <DragIndicatorIcon
            {...attributes}
            {...listeners}
            style={{
              width: "20px",
              height: "20px",
              marginRight: "8px",
              cursor: isDragging ? "grabbing" : "grab",
              color: "#999",
            }}
          />
          <span style={{ flex: 1 }}>{text}</span>
        </span>
        <IconButton style={{ padding: "0" }} onClick={() => onRemove(text)}>
          <CloseIcon style={{ width: "15px", height: "15px" }} />
        </IconButton>
      </span>
    </li>
  )
})

interface FilePathProps {
  path: string | string[]
  setSelectedFilePath: (path: string[] | string) => void
}

const FilePathSelectedListView = memo(function FilePathSelectedListView({
  path,
  setSelectedFilePath,
}: FilePathProps) {
  const sensors = useSensors(
    useSensor(PointerSensor, {
      activationConstraint: {
        distance: 8,
      },
    }),
    useSensor(KeyboardSensor, {
      coordinateGetter: sortableKeyboardCoordinates,
    }),
  )

  const handleRemoveFile = (fileToRemove: string) => {
    if (Array.isArray(path)) {
      setSelectedFilePath(path.filter((file) => file !== fileToRemove))
    }
  }

  const handleDragEnd = (event: DragEndEvent) => {
    const { active, over } = event

    if (active.id !== over?.id && Array.isArray(path)) {
      const oldIndex = path.indexOf(active.id as string)
      const newIndex = path.indexOf(over?.id as string)

      if (oldIndex !== -1 && newIndex !== -1) {
        setSelectedFilePath(arrayMove(path, oldIndex, newIndex))
      }
    }
  }

  return (
    <Typography variant="subtitle2">
      {path ? (
        Array.isArray(path) ? (
          <DndContext
            sensors={sensors}
            collisionDetection={closestCenter}
            onDragEnd={handleDragEnd}
          >
            <SortableContext
              items={path}
              strategy={verticalListSortingStrategy}
            >
              <ul style={{ padding: 0, margin: 0, listStyleType: "none" }}>
                {path.map((text) => (
                  <SortableItem
                    key={text}
                    id={text}
                    text={text}
                    onRemove={handleRemoveFile}
                  />
                ))}
              </ul>
            </SortableContext>
          </DndContext>
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
