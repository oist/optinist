import { styled } from '@mui/material'

const Input = styled('input', {
  shouldForwardProp: (props) => props !== 'error',
})<{ error?: boolean }>(({ error }) => {
  return {
    width: 250,
    height: 24,
    borderRadius: 4,
    border: '1px solid',
    borderColor: error ? 'red' : '#d9d9d9',
    padding: '5px 10px',
    marginBottom: 15,
    transition: 'all 0.3s',
    outline: 'none',
    ':focus, :hover': {
      borderColor: '#1677ff',
    },
  }
})

export default Input
