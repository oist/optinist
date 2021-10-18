import React from 'react'
import { useSelector } from 'react-redux'
import Button from '@material-ui/core/Button'
import PlayArrowIcon from '@material-ui/icons/PlayArrow'
import StopIcon from '@material-ui/icons/Stop'
import axios from 'axios'
import { FlowElement } from 'react-flow-renderer'
import { flowElementsSelector } from 'redux/slice/Element/ElementSelector'

export const ToolBar = React.memo(() => {
  const flowElements = useSelector(flowElementsSelector)

  const onRunBtnClick = () => {
    var flowList: any[] = []
    flowElements.forEach((edge: FlowElement) => {
      if ('source' in edge) {
        var node: FlowElement = flowElements.find((e) => e.id === edge.source)!
        flowList.push(node.data)
      }
    })
    // console.log(flowList)
    axios.post('http://localhost:8000/api/run', flowList).then((res) => {
      var message = res.data
      // console.log(message)
    })
    // axios.get('http://localhost:8000/api/run').then((res) => {
    //   var message = res.data
    //   console.log(message)
    // })
  }

  return (
    <>
      <Button
        className="ctrl_btn"
        variant="contained"
        color="primary"
        endIcon={<PlayArrowIcon />}
        onClick={onRunBtnClick}
      >
        run
      </Button>
      <Button
        className="ctrl_btn"
        variant="contained"
        color="secondary"
        endIcon={<StopIcon />}
      >
        stop
      </Button>
    </>
  )
})
