import MuiAccordion from '@material-ui/core/Accordion'
import { withStyles } from '@material-ui/core/styles'

export const Accordion = withStyles((theme) => ({
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
