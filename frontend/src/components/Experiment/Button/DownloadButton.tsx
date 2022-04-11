import React, { useState, useRef } from 'react'
import IconButton from '@mui/material/IconButton'
import { ExperimentUidContext } from '../ExperimentTable'
import { BASE_URL } from 'const/API'
import axios from 'axios'
import GetAppIcon from '@mui/icons-material/GetApp'

export const NWBDownloadButton = React.memo<{
  name: string
}>(({ name }) => {
  const uid = React.useContext(ExperimentUidContext)
  const ref = useRef<HTMLAnchorElement | null>(null)
  const [url, setFileUrl] = useState<string>()

  const onClick = async () => {
    try {
      const response = await axios.get(
        `${BASE_URL}/experiments/download/nwb/${uid}`,
        {
          responseType: 'blob',
        },
      )
      const url = URL.createObjectURL(new Blob([response.data]))
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
      <a href={url} download={`${name}.nwb`} className="hidden" ref={ref} />
    </>
  )
})

export const ConfigDownloadButton = React.memo(() => {
  const uid = React.useContext(ExperimentUidContext)
  const ref = useRef<HTMLAnchorElement | null>(null)
  const [url, setFileUrl] = useState<string>()

  const onClick = async () => {
    try {
      const response = await axios.get(
        `${BASE_URL}/experiments/download/config/${uid}`,
        {
          responseType: 'blob',
        },
      )
      const url = URL.createObjectURL(new Blob([response.data]))
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
      <a href={url} download={`config.yaml`} className="hidden" ref={ref} />
    </>
  )
})
