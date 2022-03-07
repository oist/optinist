import React from 'react'
import { useSelector, useDispatch } from 'react-redux'

import Button from '@mui/material/Button'
import ButtonGroup from '@mui/material/ButtonGroup'
import ArrowDropDownIcon from '@mui/icons-material/ArrowDropDown'
import CloseIcon from '@mui/icons-material/Close'
import ClickAwayListener from '@mui/material/ClickAwayListener'
import Grow from '@mui/material/Grow'
import Paper from '@mui/material/Paper'
import Popper from '@mui/material/Popper'
import MenuItem from '@mui/material/MenuItem'
import MenuList from '@mui/material/MenuList'
import TextField from '@mui/material/TextField'
import Dialog from '@mui/material/Dialog'
import DialogActions from '@mui/material/DialogActions'
import DialogContent from '@mui/material/DialogContent'
import DialogTitle from '@mui/material/DialogTitle'

import { useSnackbar } from 'notistack'

import { UseRunPipelineReturnType } from 'store/slice/Pipeline/PipelineHook'
import {
  RUN_BTN_LABELS,
  RUN_BTN_OPTIONS,
  RUN_BTN_TYPE,
} from 'store/slice/Pipeline/PipelineType'
import { selectPipelineRunBtn } from 'store/slice/Pipeline/PipelineSelectors'
import { setRunBtnOption } from 'store/slice/Pipeline/PipelineSlice'

export const RunButtons = React.memo<UseRunPipelineReturnType>((props) => {
  const {
    uid,
    isStartedSuccess,
    filePathIsUndefined,
    handleCancelPipeline,
    handleRunPipeline,
    handleRunPipelineByUid,
  } = props

  const dispatch = useDispatch()

  const runBtnOption = useSelector(selectPipelineRunBtn)

  const [dialogOpen, setDialogOpen] = React.useState(false)
  const { enqueueSnackbar } = useSnackbar()
  const handleClick = () => {
    if (!filePathIsUndefined) {
      if (runBtnOption === RUN_BTN_OPTIONS.RUN_NEW) {
        setDialogOpen(true)
      } else {
        handleRunPipelineByUid()
      }
    } else {
      enqueueSnackbar('please select input file', { variant: 'error' })
    }
  }
  const onClickDialogRun = (name: string) => {
    handleRunPipeline(name)
    setDialogOpen(false)
  }
  const onClickCancel = () => {
    handleCancelPipeline()
  }
  const [menuOpen, setMenuOpen] = React.useState(false)
  const anchorRef = React.useRef<HTMLDivElement>(null)

  const handleMenuItemClick = (
    event: React.MouseEvent<HTMLLIElement, MouseEvent>,
    option: RUN_BTN_TYPE,
  ) => {
    dispatch(setRunBtnOption({ runBtnOption: option }))
    setMenuOpen(false)
  }
  const handleToggle = () => {
    setMenuOpen((prevOpen) => !prevOpen)
  }
  const handleClose = (event: Event) => {
    if (
      anchorRef.current &&
      anchorRef.current.contains(event.target as HTMLElement)
    ) {
      return
    }
    setMenuOpen(false)
  }
  const uidExists = uid != null
  return (
    <>
      <ButtonGroup
        sx={{
          margin: (theme) => theme.spacing(1),
        }}
        variant="contained"
        ref={anchorRef}
        disabled={isStartedSuccess}
      >
        <Button onClick={handleClick}>{RUN_BTN_LABELS[runBtnOption]}</Button>
        <Button size="small" onClick={handleToggle}>
          <ArrowDropDownIcon />
        </Button>
      </ButtonGroup>
      <Popper
        open={menuOpen}
        anchorEl={anchorRef.current}
        role={undefined}
        transition
        disablePortal
      >
        {({ TransitionProps, placement }) => (
          <Grow
            {...TransitionProps}
            style={{
              transformOrigin:
                placement === 'bottom' ? 'center top' : 'center bottom',
            }}
          >
            <Paper>
              <ClickAwayListener onClickAway={handleClose}>
                <MenuList>
                  {Object.values(RUN_BTN_OPTIONS).map((option) => (
                    <MenuItem
                      key={option}
                      disabled={
                        !uidExists && option === RUN_BTN_OPTIONS.RUN_ALREADY
                      }
                      selected={option === runBtnOption}
                      onClick={(event) => handleMenuItemClick(event, option)}
                    >
                      {RUN_BTN_LABELS[option]}
                    </MenuItem>
                  ))}
                </MenuList>
              </ClickAwayListener>
            </Paper>
          </Grow>
        )}
      </Popper>
      <Button
        variant="outlined"
        endIcon={<CloseIcon />}
        onClick={onClickCancel}
        sx={{
          margin: (theme) => theme.spacing(1),
        }}
      >
        Cancel
      </Button>
      <RunDialog
        open={dialogOpen}
        handleRun={onClickDialogRun}
        handleClose={() => setDialogOpen(false)}
      />
    </>
  )
})

const RunDialog = React.memo<{
  open: boolean
  handleRun: (name: string) => void
  handleClose: () => void
}>(({ open, handleClose, handleRun }) => {
  const [name, setName] = React.useState('')
  const [error, setError] = React.useState<string | null>(null)
  const onClickRun = () => {
    if (name !== '') {
      handleRun(name)
    } else {
      setError('name is empty')
    }
  }
  const onChangeName = (event: React.ChangeEvent<HTMLInputElement>) => {
    setName(event.target.value)
    if (event.target.value !== '') {
      setError(null)
    }
  }
  return (
    <Dialog open={open} onClose={handleClose}>
      <DialogTitle>Name and run flowchart</DialogTitle>
      <DialogContent>
        <TextField
          label="name"
          autoFocus
          margin="dense"
          fullWidth
          variant="standard"
          onChange={onChangeName}
          error={error != null}
          helperText={error}
        />
      </DialogContent>
      <DialogActions>
        <Button onClick={handleClose} color="inherit" variant="outlined">
          Cancel
        </Button>
        <Button onClick={onClickRun} color="primary" variant="outlined">
          Run
        </Button>
      </DialogActions>
    </Dialog>
  )
})
