import { useDispatch, useSelector } from 'react-redux'
import Tabs, { tabsClasses } from '@mui/material/Tabs'
import Tab from '@mui/material/Tab'
import Dialog from '@mui/material/Dialog'
import DialogTitle from '@mui/material/DialogTitle'
import DialogContent from '@mui/material/DialogContent'
import IconButton from '@mui/material/IconButton'
import CloseIcon from '@mui/icons-material/Close'
import Box from '@mui/material/Box'
import React from 'react'

import { arrayEqualityFn } from 'utils/EqualityUtils'
import { selectAlgorithmName } from 'store/slice/AlgorithmNode/AlgorithmNodeSelectors'
import {
  selectPipelineNodeResultOutputKeyList,
  selectPipelineNodeResultOutputFileDataType,
  selectPipelineNodeResultOutputFilePath,
} from 'store/slice/Pipeline/PipelineSelectors'
import { selectVisualizeItemIdForWorkflowDialog } from 'store/slice/VisualizeItem/VisualizeItemSelectors'
import {
  addItemForWorkflowDialog,
  deleteAllItemForWorkflowDialog,
} from 'store/slice/VisualizeItem/VisualizeItemSlice'
import { DisplayDataItem } from 'components/Workspace/Visualize/DisplayDataItem'

export const AlgorithmOutputDialog = React.memo<{
  open: boolean
  onClose: () => void
  nodeId: string
}>(({ open, onClose, nodeId }) => {
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

const TitleWithCloseButton = React.memo<{
  onClose: () => void
  nodeId: string
}>(({ nodeId, onClose }) => {
  const nodeName = useSelector(selectAlgorithmName(nodeId))
  return (
    <DialogTitle sx={{ m: 0, p: 2 }}>
      Output of {nodeName}
      <IconButton
        onClick={onClose}
        sx={{
          position: 'absolute',
          right: 8,
          top: 10,
        }}
      >
        <CloseIcon />
      </IconButton>
    </DialogTitle>
  )
})

const OutputViewer = React.memo<{ nodeId: string }>(({ nodeId }) => {
  const outputKeyList = useSelector(
    selectPipelineNodeResultOutputKeyList(nodeId),
    arrayEqualityFn,
  )
  const [selectedOutoutKey, setSelectedOutputKey] = React.useState(
    outputKeyList[0],
  )
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

const OutputSelectTabs = React.memo<{
  selectedOutoutKey: string
  outputKeyList: string[]
  onSelectOutput: (selectedKey: string) => void
}>(({ selectedOutoutKey, outputKeyList, onSelectOutput }) => {
  const handleChange = (event: React.SyntheticEvent, newValue: string) => {
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
          '&.Mui-disabled': { opacity: 0.3 },
        },
      }}
    >
      {outputKeyList.map((outputKey) => (
        <Tab
          value={outputKey}
          label={outputKey}
          sx={{
            textTransform: 'none',
          }}
        />
      ))}
    </Tabs>
  )
})

const DisplayDataView = React.memo<{ nodeId: string; outputKey: string }>(
  ({ nodeId, outputKey }) => {
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
    React.useEffect(() => {
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
  },
)
