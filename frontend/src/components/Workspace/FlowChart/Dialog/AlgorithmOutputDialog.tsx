import { memo, SyntheticEvent, useEffect, useState } from "react"
import { useDispatch, useSelector } from "react-redux"

import CloseIcon from "@mui/icons-material/Close"
import Box from "@mui/material/Box"
import Dialog from "@mui/material/Dialog"
import DialogContent from "@mui/material/DialogContent"
import DialogTitle from "@mui/material/DialogTitle"
import IconButton from "@mui/material/IconButton"
import Tab from "@mui/material/Tab"
import Tabs, { tabsClasses } from "@mui/material/Tabs"

import { DisplayDataItem } from "components/Workspace/Visualize/DisplayDataItem"
import { selectAlgorithmName } from "store/slice/AlgorithmNode/AlgorithmNodeSelectors"
import { NodeIdProps } from "store/slice/FlowElement/FlowElementType"
import {
  selectPipelineNodeResultOutputKeyList,
  selectPipelineNodeResultOutputFileDataType,
  selectPipelineNodeResultOutputFilePath,
} from "store/slice/Pipeline/PipelineSelectors"
import { selectVisualizeItemIdForWorkflowDialog } from "store/slice/VisualizeItem/VisualizeItemSelectors"
import {
  addItemForWorkflowDialog,
  deleteAllItemForWorkflowDialog,
} from "store/slice/VisualizeItem/VisualizeItemSlice"
import { arrayEqualityFn } from "utils/EqualityUtils"

interface AlgorithmOutputDialogProps {
  open: boolean
  onClose: () => void
  nodeId: string
}

export const AlgorithmOutputDialog = memo(function AlgorithmOutputDialog({
  open,
  onClose,
  nodeId,
}: AlgorithmOutputDialogProps) {
  const dispatch = useDispatch()
  const closeFn = () => {
    onClose()
    dispatch(deleteAllItemForWorkflowDialog())
  }
  return (
    <Dialog open={open} onClose={closeFn} fullWidth>
      <TitleWithCloseButton onClose={closeFn} nodeId={nodeId} />
      <DialogContent
        dividers
        sx={{
          pt: 1,
          px: 2,
        }}
      >
        {open && <OutputViewer nodeId={nodeId} />}
      </DialogContent>
    </Dialog>
  )
})

interface TitleWithCloseButtonProps extends NodeIdProps {
  onClose: () => void
}

const TitleWithCloseButton = memo(function TitleWithCloseButtonProps({
  nodeId,
  onClose,
}: TitleWithCloseButtonProps) {
  const nodeName = useSelector(selectAlgorithmName(nodeId))
  return (
    <DialogTitle sx={{ m: 0, p: 2 }}>
      Output of {nodeName}
      <IconButton
        onClick={onClose}
        sx={{
          position: "absolute",
          right: 8,
          top: 10,
        }}
      >
        <CloseIcon />
      </IconButton>
    </DialogTitle>
  )
})

const OutputViewer = memo(function OutputViewer({ nodeId }: NodeIdProps) {
  const outputKeyList = useSelector(
    selectPipelineNodeResultOutputKeyList(nodeId),
    arrayEqualityFn,
  )
  const [selectedOutoutKey, setSelectedOutputKey] = useState(outputKeyList[0])
  return (
    <>
      <OutputSelectTabs
        outputKeyList={outputKeyList}
        selectedOutoutKey={selectedOutoutKey}
        onSelectOutput={setSelectedOutputKey}
      />
      <DisplayDataView nodeId={nodeId} outputKey={selectedOutoutKey} />
    </>
  )
})

interface OutputSelectTabsProps {
  selectedOutoutKey: string
  outputKeyList: string[]
  onSelectOutput: (selectedKey: string) => void
}

const OutputSelectTabs = memo(function OutputSelectTabs({
  selectedOutoutKey,
  outputKeyList,
  onSelectOutput,
}: OutputSelectTabsProps) {
  const handleChange = (event: SyntheticEvent, newValue: string) => {
    onSelectOutput(newValue)
  }
  return (
    <Tabs
      value={selectedOutoutKey}
      onChange={handleChange}
      variant="scrollable"
      scrollButtons="auto"
      sx={{
        [`& .${tabsClasses.scrollButtons}`]: {
          "&.Mui-disabled": { opacity: 0.3 },
        },
      }}
    >
      {outputKeyList.map((outputKey) => (
        <Tab
          key={outputKey}
          value={outputKey}
          label={outputKey}
          sx={{
            textTransform: "none",
          }}
        />
      ))}
    </Tabs>
  )
})

interface DisplayDataViewProps extends NodeIdProps {
  outputKey: string
}

const DisplayDataView = memo(function DisplayDataView({
  nodeId,
  outputKey,
}: DisplayDataViewProps) {
  const dispatch = useDispatch()
  const filePath = useSelector(
    selectPipelineNodeResultOutputFilePath(nodeId, outputKey),
  )
  const dataType = useSelector(
    selectPipelineNodeResultOutputFileDataType(nodeId, outputKey),
  )
  const itemId = useSelector(
    selectVisualizeItemIdForWorkflowDialog(nodeId, filePath, dataType),
  )
  useEffect(() => {
    if (itemId === null) {
      dispatch(addItemForWorkflowDialog({ nodeId, filePath, dataType }))
    }
  }, [dispatch, nodeId, filePath, dataType, itemId])
  return (
    <Box
      sx={{
        mx: 1,
        my: 2,
      }}
    >
      {itemId != null && <DisplayDataItem itemId={itemId} />}
    </Box>
  )
})
