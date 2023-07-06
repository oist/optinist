import { ChangeEvent, FC, FocusEvent, useState } from 'react'
import { Box, styled, Typography } from '@mui/material'
import VisibilityOffIcon from '@mui/icons-material/VisibilityOff'
import VisibilityIcon from '@mui/icons-material/Visibility'
import Input from 'components/common/Input'

const style: object = {
  position: 'absolute',
  right: 5,
  top: 8,
  fontSize: 20,
  cursor: 'pointer',
  color: 'rgba(0,0,0,0.6)',
}

type InputPasswordProps = {
  onChange?: (event: ChangeEvent<HTMLInputElement>) => void
  error?: string
  name?: string
  placeholder?: string
  onBlur?:  (event: FocusEvent<HTMLInputElement>) => void
}

const InputPassword: FC<InputPasswordProps> = ({ error, ...p }) => {
  const [type, setType] = useState('password')

  const onShow = () => {
    setType('text')
  }

  const onHidden = () => {
    setType('password')
  }

  return (
    <Box sx={{ position: 'relative' }}>
      <Input {...p} type={type} />
      {type === 'password' ? (
        <VisibilityIcon style={style} onClick={onShow} />
      ) : (
        <VisibilityOffIcon style={style} onClick={onHidden} />
      )}
      <TextError>{error}</TextError>
    </Box>
  )
}

const TextError = styled(Typography)({
  fontSize: 12,
  minHeight: 18,
  color: 'red',
  lineHeight: '14px',
  marginTop: -14,
  wordBreak: 'break-word',
})

export default InputPassword
