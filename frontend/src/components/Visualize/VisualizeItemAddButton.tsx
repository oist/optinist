import React from 'react'
import { useDispatch } from 'react-redux'
import { styled } from '@mui/material/styles'
import Box from '@mui/material/Box'
import Paper from '@mui/material/Paper'
import AddIcon from '@mui/icons-material/Add'
import Button from '@mui/material/Button'
import { addInitialItem } from 'store/slice/VisualizeItem/VisualizeItemSlice'

export const VisualizeItemAddButton: React.FC = () => {
  const dispatch = useDispatch()
  const onClick = () => {
    dispatch(addInitialItem())
  }
  return (
    <StyledPaper elevation={1} variant="outlined">
      <Box
        display="flex"
        justifyContent="center"
        alignItems="center"
        height="100%"
      >
        <StyledButton onClick={onClick}>
          <AddIcon fontSize="large" color="primary" />
        </StyledButton>
      </Box>
    </StyledPaper>
  )
}

const StyledPaper = styled(Paper)(({ theme }) => ({
  width: 260,
  height: 255,
  border: 'dashed',
  borderWidth: 2,
  borderColor: theme.palette.divider,
  margin: theme.spacing(1),
}))

const StyledButton = styled(Button)({
  width: '100%',
  height: '100%',
})
