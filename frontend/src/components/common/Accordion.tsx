import MuiAccordion, { AccordionProps } from '@mui/material/Accordion'
import { styled } from '@mui/material/styles'

export const Accordion = styled((props: AccordionProps) => (
  <MuiAccordion disableGutters elevation={0} square {...props} />
))(({ theme }) => ({
  border: '1px solid',
  borderColor: theme.palette.divider,
  boxShadow: 'none',
  '&:not(:last-child)': {
    borderBottom: 0,
  },
  '&:before': {
    display: 'none',
  },
}))
