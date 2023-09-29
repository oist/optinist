import { InputProps, styled, Typography } from '@mui/material'

interface InputErrorProps extends InputProps {
  errorMessage: string
  value?: string
}

const InputError = ({
  errorMessage,
  onChange,
  value,
  type,
  onBlur,
  name,
}: InputErrorProps) => {
  return (
    <>
      <Input
        autoComplete="Off"
        error={!!errorMessage}
        onChange={onChange}
        value={value}
        type={type}
        name={name}
        onBlur={onBlur}
      />
      <TextError>{errorMessage}</TextError>
    </>
  )
}

const Input = styled('input', {
  shouldForwardProp: (props) => props !== 'error',
})<{ error: boolean }>(({ error }) => {
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

const TextError = styled(Typography)({
  fontSize: 12,
  minHeight: 18,
  color: 'red',
  lineHeight: '14px',
  margin: '-14px 0px 0px 305px',
  wordBreak: 'break-word',
})

export default InputError
