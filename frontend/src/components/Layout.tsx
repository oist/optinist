import React from 'react'
import AppBar from '@material-ui/core/AppBar'
import Tabs from '@material-ui/core/Tabs'
import Tab from '@material-ui/core/Tab'
import Typography from '@material-ui/core/Typography'
import { makeStyles } from '@material-ui/core/styles'
import Toolbar from '@material-ui/core/Toolbar'
import FlowChart from './FlowChart/FlowChart'
import Visualize from './Visualize/Visualize'
import { useLazyRunPipelineQuery } from 'api/Run/Run'
import { RunPipeLineContext } from './Visualize/DataContext'

const Layout: React.FC = () => {
  const classes = useStyles()
  const [value, setValue] = React.useState(0)
  const handleChange = (event: React.ChangeEvent<{}>, newValue: number) => {
    setValue(newValue)
  }
  const [runPipeLine, result] = useLazyRunPipelineQuery()

  return (
    <div className={classes.root}>
      <AppBar position="fixed" className={classes.appBar}>
        <Toolbar variant="dense">
          <Typography variant="h6">OPTINIST</Typography>
          <Tabs
            className={classes.tab}
            value={value}
            onChange={handleChange}
            centered
          >
            <Tab label="Flow Chart" {...a11yProps(0)} />
            <Tab label="Visualize" {...a11yProps(1)} />
          </Tabs>
        </Toolbar>
      </AppBar>
      <TabPanel value={value} index={0}>
        <RunPipeLineContext.Provider value={{ runPipeLine, result }}>
          <FlowChart />
        </RunPipeLineContext.Provider>
      </TabPanel>
      <TabPanel value={value} index={1}>
        <Visualize />
      </TabPanel>
    </div>
  )
}

const useStyles = makeStyles((theme) => ({
  root: {
    flexGrow: 1,
    backgroundColor: theme.palette.background.paper,
  },
  appBar: {
    zIndex: theme.zIndex.drawer + 1,
  },
  tab: {
    width: '100%',
  },
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
