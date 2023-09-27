import {
  MenuItem,
  Select,
  SelectChangeEvent,
  styled,
  Typography,
} from '@mui/material'
import { FC, FocusEvent } from 'react'

type SelectErrorProps = {
  value?: string
  onChange?: (value: SelectChangeEvent, child: React.ReactNode) => void
  onBlur?: (event: FocusEvent<HTMLInputElement | HTMLTextAreaElement>) => void
  errorMessage: string
  name?: string
  options: string[]
}

const SelectError: FC<SelectErrorProps> = ({
  value,
  onChange,
  onBlur,
  errorMessage,
  options,
  name,
}) => {
  return (
    <>
      <SelectModal
        name={name}
        value={value}
        onChange={
          onChange as (
            value: SelectChangeEvent<unknown>,
            child: React.ReactNode,
          ) => void
        }
        onBlur={onBlur}
        error={!!errorMessage}
      >
        {options.map((item: string) => {
          return (
            <MenuItem key={item} value={item}>
              {item}
            </MenuItem>
          )
        })}
      </SelectModal>
      <TextError>{errorMessage}</TextError>
    </>
  )
}

const SelectModal = styled(Select, {
  shouldForwardProp: (props) => props !== 'error',
})<{ error: boolean }>(({ theme, error }) => ({
  width: 272,
  marginBottom: '15px',
  border: '1px solid #d9d9d9',
  borderColor: error ? 'red' : '#d9d9d9',
  borderRadius: 4,
}))

const TextError = styled(Typography)({
  fontSize: 12,
  minHeight: 18,
  color: 'red',
  lineHeight: '14px',
  margin: '-14px 0px 0px 305px',
  wordBreak: 'break-word',
})
export default SelectError
