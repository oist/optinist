import React from 'react'
import AppBar from '@mui/material/AppBar'
import Tabs from '@mui/material/Tabs'
import Tab from '@mui/material/Tab'
import Typography from '@mui/material/Typography'
import Toolbar from '@mui/material/Toolbar'
import { styled } from '@mui/material/styles'
import FlowChart from './FlowChart/FlowChart'
import Visualize from './Visualize/Visualize'
import Experiment from './Experiment/Experiment'
import { useRunPipeline } from 'store/slice/Pipeline/PipelineHook'

const Layout: React.FC = () => {
  const [value, setValue] = React.useState(0)
  const handleChange = (event: React.ChangeEvent<{}>, newValue: number) => {
    setValue(newValue)
  }

  const runPipeline = useRunPipeline() // タブ切り替えによって結果取得処理が止まってしまうのを回避するため、タブの親レイヤーで呼び出している

  return (
    <RootDiv>
      <StyledAppBar position="fixed">
        <Toolbar variant="dense">
          <Typography variant="h6">OPTINIST</Typography>
          <Tabs
            sx={{ width: '100%' }}
            value={value}
            onChange={handleChange}
            centered
            textColor="inherit"
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
}))

const StyledAppBar = styled(AppBar)(({ theme }) => ({
  zIndex: theme.zIndex.drawer + 1,
}))

interface TabPanelProps {
  children?: React.ReactNode
  index: any
  value: any
}

function TabPanel(props: TabPanelProps) {
  const { children, value, index, ...other } = props

  return (
    <div
      role="tabpanel"
      hidden={value !== index}
      id={`simple-tabpanel-${index}`}
      aria-labelledby={`simple-tab-${index}`}
      {...other}
    >
      {value === index && <Typography>{children}</Typography>}
    </div>
  )
}

function a11yProps(index: any) {
  return {
    id: `simple-tab-${index}`,
    'aria-controls': `simple-tabpanel-${index}`,
  }
}

export default Layout
