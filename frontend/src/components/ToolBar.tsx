import react from 'react'
import Button from '@material-ui/core/Button'
import PlayArrowIcon from '@material-ui/icons/PlayArrow'
import StopIcon from '@material-ui/icons/Stop'
import axios from 'axios'

const ToolBar: react.FC = () => {
  const onRunBtnClick = (event: React.MouseEvent<Element, MouseEvent>) => {
    console.log(event)
    console.log('run')
    axios.get('http://localhost:8000/run').then((res) => {
      var message = res.data
      console.log(message)
    })
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
