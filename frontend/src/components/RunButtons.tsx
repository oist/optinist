import React from 'react'

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

const OPTIONS = {
  RUN_NEW: 1,
  RUN_ALREADY: 2,
} as const

export type OPTION_TYPE = typeof OPTIONS[keyof typeof OPTIONS]

const OPTIONS_LABELS = {
  [OPTIONS.RUN_NEW]: 'RUN(NEW)',
  [OPTIONS.RUN_ALREADY]: 'RUN(ALREADY)',
} as const

export const RunButtons = React.memo<UseRunPipelineReturnType>((props) => {
  const {
    uid,
    isStartedSuccess,
    filePathIsUndefined,
    handleCancelPipeline,
    handleRunPipeline,
    handleRunPipelineByUid,
  } = props

  const [dialogOpen, setDialogOpen] = React.useState(false)
  const { enqueueSnackbar } = useSnackbar()
  const handleClick = () => {
    if (!filePathIsUndefined) {
      if (selectedOption === OPTIONS.RUN_NEW) {
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
  const [selectedOption, setSelectedOption] = React.useState<OPTION_TYPE>(
    OPTIONS.RUN_NEW,
  )
  const handleMenuItemClick = (
    event: React.MouseEvent<HTMLLIElement, MouseEvent>,
    option: OPTION_TYPE,
  ) => {
    setSelectedOption(option)
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
        <Button onClick={handleClick}>{OPTIONS_LABELS[selectedOption]}</Button>
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
                  {Object.values(OPTIONS).map((option) => (
                    <MenuItem
                      key={option}
                      disabled={!uidExists && option === OPTIONS.RUN_ALREADY}
                      selected={option === selectedOption}
                      onClick={(event) => handleMenuItemClick(event, option)}
                    >
                      {OPTIONS_LABELS[option]}
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
