import { FC, SyntheticEvent, useCallback } from "react"
import { useDispatch, useSelector } from "react-redux"
import { useNavigate } from "react-router-dom"

import ExitToAppIcon from "@mui/icons-material/ExitToApp"
import Tab from "@mui/material/Tab"
import Tabs from "@mui/material/Tabs"

import { StatusROI } from "components/Workspace/Visualize/Plot/ImagePlot"
import { cancelRoi } from "store/slice/DisplayData/DisplayDataActions"
import { selectPipelineIsStartedSuccess } from "store/slice/Pipeline/PipelineSelectors"
import { selectModeStandalone } from "store/slice/Standalone/StandaloneSeclector"
import {
  selectActiveTab,
  selectCurrentWorkspaceId,
  selectRoiFilePathCancel,
  selectStatusRoiCancel,
} from "store/slice/Workspace/WorkspaceSelector"
import { setActiveTab } from "store/slice/Workspace/WorkspaceSlice"
import { AppDispatch } from "store/store"

const WorkspaceTabs: FC = () => {
  const dispatch = useDispatch<AppDispatch>()
  const activeTab = useSelector(selectActiveTab)
  const roiFilePath = useSelector(selectRoiFilePathCancel)
  const workspaceId = useSelector(selectCurrentWorkspaceId)
  const statusRoi = useSelector(selectStatusRoiCancel)
  const navigate = useNavigate()
  const isStandalone = useSelector(selectModeStandalone)

  const handleChange = useCallback(
    (event: SyntheticEvent, newValue: number) => {
      dispatch(setActiveTab(newValue))
      if (!statusRoi || !roiFilePath || !workspaceId) return
      const checkCancel = !Object.keys(statusRoi).every(
        (key) => statusRoi[key as keyof StatusROI].length === 0,
      )
      if (newValue !== 1 && checkCancel) {
        dispatch(cancelRoi({ path: roiFilePath, workspaceId }))
      }
    },
    //eslint-disable-next-line
    [roiFilePath, workspaceId, statusRoi],
  )

  const isPending = useSelector(selectPipelineIsStartedSuccess)

  return (
    <>
      <Tabs
        sx={{ width: "100%" }}
        value={activeTab}
        onChange={handleChange}
        centered
        textColor="primary"
      >
        <Tab label="Workflow" {...a11yProps(0)} />
        <Tab disabled={!!isPending} label="Visualize" {...a11yProps(1)} />
        <Tab label="Record" {...a11yProps(2)} />
        {!isStandalone && (
          <Tab
            label={
              <span>
                <span style={{ display: "flex", alignItems: "center" }}>
                  <ExitToAppIcon
                    fontSize="small"
                    style={{ marginRight: "0.2rem" }}
                  />
                  <span style={{ fontSize: "smaller" }}>Workspaces</span>
                </span>
              </span>
            }
            onClick={() => navigate("/console/workspaces")}
          />
        )}
      </Tabs>
    </>
  )
}

function a11yProps(index: number | string) {
  return {
    id: `simple-tab-${index}`,
    "aria-controls": `simple-tabpanel-${index}`,
  }
}

export default WorkspaceTabs
