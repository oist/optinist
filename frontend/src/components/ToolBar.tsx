import Button from '@material-ui/core/Button'
import PlayArrowIcon from '@material-ui/icons/PlayArrow'
import StopIcon from '@material-ui/icons/Stop'
import react from 'react'

const ToolBar: react.FC = () => {
  return (
    <>
      <Button
        className="ctrl_btn"
        variant="contained"
        color="primary"
        endIcon={<PlayArrowIcon />}
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
