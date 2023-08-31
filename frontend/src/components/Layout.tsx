import React from 'react'
import AppBar from '@mui/material/AppBar'
import Tabs from '@mui/material/Tabs'
import Tab from '@mui/material/Tab'
import Typography from '@mui/material/Typography'
import Toolbar from '@mui/material/Toolbar'
import IconButton from '@mui/material/IconButton'
import Tooltip from '@mui/material/Tooltip'
import { styled } from '@mui/material/styles'
import GitHubIcon from '@mui/icons-material/GitHub'
import MenuBookIcon from '@mui/icons-material/MenuBook'

import { useRunPipeline } from 'store/slice/Pipeline/PipelineHook'
import FlowChart from './FlowChart/FlowChart'
import Visualize from './Visualize/Visualize'
import Experiment from './Experiment/Experiment'
import optinistLogo from './optinist.png'

const Layout: React.FC = () => {
  const [value, setValue] = React.useState(0)
  const handleChange = (event: React.ChangeEvent<{}>, newValue: number) => {
    setValue(newValue)
  }

  const runPipeline = useRunPipeline() // タブ切り替えによって結果取得処理が止まってしまうのを回避するため、タブの親レイヤーで呼び出している

  return (
    <RootDiv>
      <StyledAppBar position="fixed" color="inherit">
        <Toolbar variant="dense">
          <img src={optinistLogo} alt="optinist" width={75} height={50} />
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
          <Tooltip title="GitHub repository">
            <IconButton
              sx={{ mr: 1 }}
              href="https://github.com/oist/optinist"
              target="_blank"
            >
              <GitHubIcon />
            </IconButton>
          </Tooltip>
          <Tooltip title="Documentation">
            <IconButton
              href="https://optinist.readthedocs.io/en/latest/"
              target="_blank"
            >
              <MenuBookIcon />
            </IconButton>
          </Tooltip>
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
  backgroundColor: '#E1DEDB',
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
      {value === index && <div>{children}</div>}
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
