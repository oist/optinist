import React, { useState, useRef } from 'react'
import IconButton from '@mui/material/IconButton'
import GetAppIcon from '@mui/icons-material/GetApp'

import {
  downloadExperimentConfigApi,
  downloadExperimentNwbApi,
} from 'api/experiments/Experiments'

import { ExperimentUidContext } from '../ExperimentTable'

export const NWBDownloadButton = React.memo<{
  name: string
}>(({ name }) => {
  const uid = React.useContext(ExperimentUidContext)
  const ref = useRef<HTMLAnchorElement | null>(null)
  const [url, setFileUrl] = useState<string>()

  const onClick = async () => {
    try {
      const responseData = await downloadExperimentNwbApi(uid)
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
      <a href={url} download={`${name}.nwb`} className="hidden" ref={ref}>
        {/* 警告が出るので空文字を入れておく */}{' '}
      </a>
    </>
  )
})

export const ConfigDownloadButton = React.memo(() => {
  const uid = React.useContext(ExperimentUidContext)
  const ref = useRef<HTMLAnchorElement | null>(null)
  const [url, setFileUrl] = useState<string>()

  const onClick = async () => {
    try {
      const responseData = await downloadExperimentConfigApi(uid)
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
