import React, { useEffect } from 'react'
import { useSelector, useDispatch } from 'react-redux'
import IconButton from '@material-ui/core/IconButton'
import CloseIcon from '@material-ui/icons/Close'
import TextField from '@material-ui/core/TextField'
import Box from '@material-ui/core/Box'
import Button from '@material-ui/core/Button'
import ExpandMoreIcon from '@material-ui/icons/ExpandMore'
import MuiAccordion from '@material-ui/core/Accordion'
import AccordionDetails from '@material-ui/core/AccordionDetails'
import AccordionSummary from '@material-ui/core/AccordionSummary'
import Typography from '@material-ui/core/Typography'
import Dialog from '@material-ui/core/Dialog'
import DialogContent from '@material-ui/core/DialogContent'
import MuiDialogTitle from '@material-ui/core/DialogTitle'
import {
  createStyles,
  Theme,
  withStyles,
  WithStyles,
} from '@material-ui/core/styles'
import { nwbListSelector } from 'store/slice/NWB/NWBSelector'
import { NWBChild, NWBNodeType } from 'store/slice/NWB/NWBType'
import { updateParam } from 'store/slice/NWB/NWB'
import { getNWBParams } from 'store/slice/NWB/NWBAction'

export const NWB = React.memo(() => {
  const [open, setOpen] = React.useState(false)
  const anchorElRef = React.useRef<HTMLButtonElement | null>(null)

  const handleClick = () => {
    setOpen((prevOpen) => !prevOpen)
  }

  const handleClose = () => {
    setOpen(false)
  }

  return (
    <div>
      <Button
        variant="contained"
        color="default"
        onClick={handleClick}
        ref={anchorElRef}
      >
        NWB setting
      </Button>
      <Dialog open={open} onClose={handleClose}>
        <DialogTitle onClose={handleClose}>NWB Parameter</DialogTitle>
        <DialogContent dividers>
          <ParamBar />
        </DialogContent>
      </Dialog>
    </div>
  )
})

const ParamBar = React.memo(() => {
  const nwbList = useSelector(nwbListSelector)

  const dispatch = useDispatch()

  useEffect(() => {
    if (Object.keys(nwbList).length === 0) {
      dispatch(getNWBParams())
    }
  })
  return (
    <div className="nwbParam">
      {Object.entries(nwbList).map(([name, node], i) => (
        <ParamComponent key={name} name={name} node={node} />
      ))}
    </div>
  )
})

const ParamComponent = React.memo<{
  name: string
  node: NWBNodeType
}>(({ name, node }) => {
  if (node.type === 'child') {
    return <ParamItem name={name} node={node} />
  } else {
    return (
      <Accordion square>
        <AccordionSummary
          expandIcon={<ExpandMoreIcon />}
          aria-controls="panel1bh-content"
          id="panel1bh-header"
        >
          {name}
        </AccordionSummary>
        <AccordionDetails>
          <div>
            {Object.entries(node.children).map(([name, childNode], i) => (
              <ParamComponent name={name} node={childNode} />
            ))}
          </div>
        </AccordionDetails>
      </Accordion>
    )
  }
})

const ParamItem = React.memo<{
  name: string
  node: NWBChild
}>(({ name, node: { value, path } }) => {
  const dispatch = useDispatch()
  const onChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    dispatch(updateParam({ paramPath: path, newValue: event.target.value }))
  }
  return (
    <Box
      sx={{
        display: 'flex',
        marginTop: 16,
        marginBottom: 16,
        alignItems: 'center',
      }}
    >
      <Box
        style={{ verticalAlign: 'middle' }}
        sx={{
          flexGrow: 1,
          width: '50%',
        }}
      >
        <Typography>{name}</Typography>
      </Box>
      <Box sx={{ width: '50%' }}>
        <TextField value={String(value)} onChange={onChange} />
      </Box>
    </Box>
  )
})

const styles = (theme: Theme) =>
  createStyles({
    root: {
      margin: 0,
      padding: theme.spacing(2),
    },
    closeButton: {
      position: 'absolute',
      right: theme.spacing(1),
      top: theme.spacing(1),
      color: theme.palette.grey[500],
    },
  })

interface DialogTitleProps extends WithStyles<typeof styles> {
  id?: string
  children: React.ReactNode
  onClose: () => void
}

const DialogTitle = withStyles(styles)((props: DialogTitleProps) => {
  const { children, classes, onClose, ...other } = props
  return (
    <MuiDialogTitle disableTypography className={classes.root} {...other}>
      <Typography variant="h6">{children}</Typography>
      {onClose ? (
        <IconButton
          aria-label="close"
          className={classes.closeButton}
          onClick={onClose}
        >
          <CloseIcon />
        </IconButton>
      ) : null}
    </MuiDialogTitle>
  )
})

const Accordion = withStyles((theme) => ({
  root: {
    border: '1px solid',
    borderColor: theme.palette.divider,
    boxShadow: 'none',
    '&:not(:last-child)': {
      borderBottom: 0,
    },
    '&:before': {
      display: 'none',
    },
    '&$expanded': {
      margin: 'auto',
    },
  },
  expanded: {},
}))(MuiAccordion)
