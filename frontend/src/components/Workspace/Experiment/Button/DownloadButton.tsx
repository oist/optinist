import React, { useState, useRef } from 'react'
import IconButton from '@mui/material/IconButton'
import GetAppIcon from '@mui/icons-material/GetApp'
import { useSelector } from 'react-redux'
import {
  downloadExperimentConfigApi,
  downloadExperimentNwbApi,
} from 'api/experiments/Experiments'

import { ExperimentUidContext } from '../ExperimentTable'
import { selectCurrentWorkspaceId } from 'store/slice/Workspace/WorkspaceSelector'

export const NWBDownloadButton = React.memo<{
  name: string
  nodeId?: string
  hasNWB: boolean
}>(({ name, nodeId, hasNWB }) => {
  const workspaceId = useSelector(selectCurrentWorkspaceId)
  const uid = React.useContext(ExperimentUidContext)
  const ref = useRef<HTMLAnchorElement | null>(null)
  const [url, setFileUrl] = useState<string>()

  const onClick = async () => {
    try {
      const responseData = await downloadExperimentNwbApi(workspaceId!, uid, nodeId)
      const url = URL.createObjectURL(new Blob([responseData]))
      setFileUrl(url)
      ref.current?.click()
      URL.revokeObjectURL(url)
    } catch (error) {
      throw new Error('Download Error')
    }
  }

  return (
    <>
      <IconButton onClick={onClick} color="primary" disabled={!hasNWB}>
        <GetAppIcon />
      </IconButton>
      <a href={url} download={`${name}.nwb`} className="hidden" ref={ref}>
        {/* 警告が出るので空文字を入れておく */}{' '}
      </a>
    </>
  )
})

export const ConfigDownloadButton = React.memo(() => {
  const workspaceId = useSelector(selectCurrentWorkspaceId)
  const uid = React.useContext(ExperimentUidContext)
  const ref = useRef<HTMLAnchorElement | null>(null)
  const [url, setFileUrl] = useState<string>()

  const onClick = async () => {
    try {
      const responseData = await downloadExperimentConfigApi(workspaceId!, uid)
      const url = URL.createObjectURL(new Blob([responseData]))
      setFileUrl(url)
      ref.current?.click()
      URL.revokeObjectURL(url)
    } catch (error) {
      throw new Error('Download Error')
    }
  }

  return (
    <>
      <IconButton onClick={onClick}>
        <GetAppIcon color="primary" />
      </IconButton>
      <a href={url} download={`config.yaml`} className="hidden" ref={ref}>
        {/* 警告が出るので空文字を入れておく */}{' '}
      </a>
    </>
  )
})
