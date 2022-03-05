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

import { UseRunPipelineReturnType } from 'store/slice/Pipeline/PipelineHook'
import { useSnackbar } from 'notistack'

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

  // タブ移動による再レンダリングするたびにスナックバーが実行されてしまう挙動を回避するために前回の値を保持
  // const [prevStatus, setPrevStatus] = React.useState(status)
  // React.useEffect(() => {
  //   if (prevStatus !== status) {
  //     if (status === RUN_STATUS.FINISHED) {
  //       enqueueSnackbar('Finished', { variant: 'success' })
  //     } else if (status === RUN_STATUS.ABORTED) {
  //       enqueueSnackbar('Aborted', { variant: 'error' })
  //     } else if (status === RUN_STATUS.CANCELED) {
  //       enqueueSnackbar('Canceled', { variant: 'info' })
  //     }
  //     setPrevStatus(status)
  //   }
  // }, [status, prevStatus, enqueueSnackbar])
  const [menuOpen, setMenuOpen] = React.useState(false)
  const anchorRef = React.useRef<HTMLDivElement>(null)
  const [selectedOption, setSelectedOption] = React.useState<OPTION_TYPE>(
    OPTIONS.RUN_NEW,
  )
  const { enqueueSnackbar } = useSnackbar()
  const handleClick = () => {
    if (!filePathIsUndefined) {
      if (selectedOption === OPTIONS.RUN_NEW) {
        handleRunPipeline()
      } else {
        handleRunPipelineByUid()
      }
    } else {
      enqueueSnackbar('please select input file', { variant: 'error' })
    }
  }
  const onClickCancel = () => {
    handleCancelPipeline()
  }
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
    </>
  )
})
