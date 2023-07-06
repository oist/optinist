import React, { useEffect } from 'react'
import { useParams } from 'react-router-dom'
import Tabs from '@mui/material/Tabs'
import Tab from '@mui/material/Tab'
import Toolbar from '@mui/material/Toolbar'
import { styled } from '@mui/material/styles'
import { useRunPipeline } from 'store/slice/Pipeline/PipelineHook'
import FlowChart from '../../components/Workspace/FlowChart/FlowChart'
import Experiment from '../../components/Workspace/Experiment/Experiment'
import { Box } from '@mui/material'
import Visualize from "../../components/Workspace/Visualize/Visualize";
import { useDispatch } from 'react-redux'
import { clearCurrentWorkspace, setCurrentWorkspace } from 'store/slice/Workspace/WorkspaceSlice'

const Workspace: React.FC = () => {
  const dispatch = useDispatch()
  const [value, setValue] = React.useState(0)
  const handleChange = (event: React.ChangeEvent<{}>, newValue: number) => {
    setValue(newValue)
  }
  const runPipeline = useRunPipeline() // タブ切り替えによって結果取得処理が止まってしまうのを回避するため、タブの親レイヤーで呼び出している

  const { workspaceId } = useParams<{ workspaceId: string }>()
  useEffect(() => {
    workspaceId && dispatch(setCurrentWorkspace(workspaceId))
    return () => {
      dispatch(clearCurrentWorkspace())
    }
  }, [workspaceId, dispatch])

  return (
    <RootDiv>
      <StyledAppBar color="inherit">
        <Toolbar variant="dense">
          <Tabs
            sx={{ width: '100%' }}
            value={value}
            onChange={handleChange}
            centered
            textColor="primary"
          >
            <Tab label="Workflow" {...a11yProps(0)} />
            <Tab label="Visualize" {...a11yProps(1)} />
            <Tab label="Record" {...a11yProps(2)} />
          </Tabs>
        </Toolbar>
      </StyledAppBar>
      <TabPanel value={value} index={0}>
        <FlowChart {...runPipeline} />
      </TabPanel>
      <TabPanel value={value} index={1}>
        <Visualize />
      </TabPanel>
      <TabPanel value={value} index={2}>
        <Experiment />
      </TabPanel>
    </RootDiv>
  )
}

const RootDiv = styled('div')(({ theme }) => ({
  flexGrow: 1,
  backgroundColor: theme.palette.background.paper,
  height: '100%',
}))

const StyledAppBar = styled(Box)(({ theme }) => ({
  zIndex: theme.zIndex.drawer + 1,
  backgroundColor: '#E1DEDB',
}))

interface TabPanelProps {
  children?: React.ReactNode
  index: number
  value: number
}

function TabPanel(props: TabPanelProps) {
  const { children, value, index, ...other } = props

  return (
    <div
      style={{ height: 'calc(100% - 58px)' }}
      role="tabpanel"
      hidden={value !== index}
      id={`simple-tabpanel-${index}`}
      aria-labelledby={`simple-tab-${index}`}
      {...other}
    >
      {value === index && <Box sx={{ height: '100%' }}>{children}</Box>}
    </div>
  )
}

function a11yProps(index: number | string) {
  return {
    id: `simple-tab-${index}`,
    'aria-controls': `simple-tabpanel-${index}`,
  }
}

export default Workspace
