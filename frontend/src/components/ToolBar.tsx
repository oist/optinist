import react from 'react'
import Button from '@material-ui/core/Button'
import PlayArrowIcon from '@material-ui/icons/PlayArrow'
import StopIcon from '@material-ui/icons/Stop'
import axios from 'axios'
import { useSelector } from 'react-redux'
import { flowElementsSelector } from 'redux/slice/Element/ElementSelector'
import { FlowElement } from 'react-flow-renderer'

const ToolBar: react.FC = () => {
  const flowElements = useSelector(flowElementsSelector)

  const onRunBtnClick = () => {
    var flowList: any[] = []
    flowElements.forEach((edge: FlowElement) => {
      if ('source' in edge) {
        var node: FlowElement = flowElements.find((e) => e.id == edge.source)!
        flowList.push(node.data)
      }
    })
    console.log(flowList)
    axios.post('http://localhost:8000/run', flowList).then((res) => {
      var message = res.data
      console.log(message)
    })
    // axios.get('http://localhost:8000/run').then((res) => {
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
}

export default ToolBar
