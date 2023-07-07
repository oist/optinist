import { FC } from 'react'
import { useDispatch, useSelector } from 'react-redux'
import Tabs from '@mui/material/Tabs'
import Tab from '@mui/material/Tab'
import { selectActiveTab } from 'store/slice/Workspace/WorkspaceSelector'
import { setActiveTab } from 'store/slice/Workspace/WorkspaceSlice'

const WorkspaceTabs: FC = () => {
  const dispatch = useDispatch()
  const activeTab = useSelector(selectActiveTab)
  const handleChange = (event: React.ChangeEvent<{}>, newValue: number) => {
    dispatch(setActiveTab(newValue))
  }

  return (
    <Tabs
      sx={{ width: '100%' }}
      value={activeTab}
      onChange={handleChange}
      centered
      textColor="primary"
    >
      <Tab label="Workflow" {...a11yProps(0)} />
      <Tab label="Visualize" {...a11yProps(1)} />
      <Tab label="Record" {...a11yProps(2)} />
    </Tabs>
  )
}

function a11yProps(index: number | string) {
  return {
    id: `simple-tab-${index}`,
    'aria-controls': `simple-tabpanel-${index}`,
  }
}

export default WorkspaceTabs
