import { Box, Button, Modal, styled, Typography } from '@mui/material'
import { ChangeEvent, FC, useState } from 'react'
import { regexIgnoreS, regexPassword } from 'const/Auth'
import InputPassword from 'components/Account/InputPassword'

type ChangePasswordModalProps = {
  onClose: () => void
  open: boolean
  onSubmit: (oldPass: string, newPass: string) => void
}

const ChangePasswordModal: FC<ChangePasswordModalProps> = ({
  onClose,
  open,
  onSubmit,
}) => {
  const [errors, setErrors] = useState<{ [key: string]: string }>({})
  const [values, setValues] = useState<{ [key: string]: string }>({})
  const onChangeValue = (
    event: ChangeEvent<HTMLInputElement>,
    validate?: Function,
  ) => {
    const { name, value } = event.target
    setValues({ ...values, [name]: value })
    if(name === 'new_password' && values.confirm_password) {
      if(!validate?.(value)) {
        setErrors({ ...errors, [name]: validate?.(value), confirm_password: value !== values.confirm_password ? 'Passwords do not match' : ''})
        return
      }
        setErrors({ ...errors, [name]: validate?.(value), confirm_password: ''})
        return
    }
    setErrors({ ...errors, [name]: validate?.(value) })
  }

  const validatePassword = (value: string): string => {
    if (!value) return 'This field is required'
    if(value.length > 255) return 'The text may not be longer than 255 characters'
    if (!regexPassword.test(value)) {
      return 'Your password must be at least 6 characters long and must contain at least one letter, number, and special character'
    }
    if(regexIgnoreS.test(value)) {
      return 'Allowed special characters (!#$%&()*+,-./@_|)'
    }
    return ''
  }

  const validateReEnter = (value: string): string => {
    if (!value) return 'This field is required'
    if (value !== values.new_password) {
      return 'Passwords do not match'
    }
    return ''
  }

  const validateReEnterWhenInputPassword = () => {
    const { reEnter, new_password } = values
    if (!new_password) setErrors((pre) => ({ ...pre, new_password: 'This field is required' }))
    if (reEnter && reEnter !== new_password) {
      setErrors((pre) => ({ ...pre, reEnter: 'Passwords do not match' }))
    }
  }

  const validateForm = () => {
    const errorNewPass = validatePassword(values.new_password)
    const errorConfirmPass = validatePassword(values.confirm_password)
    return {
      new_password: errorNewPass,
      confirm_password: errorConfirmPass
    }
  }

  const onChangePass = () => {
    const newErrors: { [key: string]: string } = validateForm()
    if(errors.new_password || errors.confirm_password) return
    if(newErrors.new_password || newErrors.confirm_password) {
      setErrors(newErrors)
      return
    }
    onSubmit(values.password, values.new_password)
  }

  const onCloseModal = () => {
    setErrors({})
    setValues({})
    onClose()
  }

  return (
    <Modal
      open={open}
      onClose={onCloseModal}
      aria-labelledby="modal-modal-title"
      aria-describedby="modal-modal-description"
    >
      <ContentDelete>
        <BoxTitle>
          <Typography sx={{ fontWeight: 600, fontSize: 18 }}>
            Change Password
          </Typography>
          <Typography style={{ fontSize: 13 }}>
            <span style={{ color: 'red' }}>*</span> is required
          </Typography>
        </BoxTitle>
        <BoxConfirm>
          <FormInline>
            <Label>
              Old Password <span style={{ color: 'red' }}>*</span>
            </Label>
            <InputPassword
              onChange={(e) => onChangeValue(e)}
              name="password"
              error={errors.password}
              onBlur={(e) => onChangeValue(e)}
              placeholder="Old Password"
            />
          </FormInline>
          <FormInline>
            <Label>
              New Password <span style={{ color: 'red' }}>*</span>
            </Label>
            <InputPassword
              onChange={(e) => onChangeValue(e, validatePassword)}
              name="new_password"
              error={errors.new_password}
              placeholder="New Password"
              onBlur={validateReEnterWhenInputPassword}
            />
          </FormInline>
          <FormInline>
            <Label>
              Confirm Password <span style={{ color: 'red' }}>*</span>
            </Label>
            <InputPassword
              onChange={(e) => onChangeValue(e, validateReEnter)}
              name="confirm_password"
              error={errors.confirm_password}
              placeholder="Confirm Password"
              onBlur={(e) => onChangeValue(e, validateReEnter)}
            />
          </FormInline>
          <Box sx={{ display: 'flex', justifyContent: 'flex-end' }}>
            <ButtonConfirm onClick={() => onChangePass()}>UPDATE</ButtonConfirm>
          </Box>
        </BoxConfirm>
        <Button onClick={onCloseModal}>
          <Typography
            sx={{
              textDecoration: 'underline',
              textTransform: 'none',
              lineHeight: '17px',
            }}
          >
            Close
          </Typography>
        </Button>
      </ContentDelete>
    </Modal>
  )
}

const BoxTitle = styled(Box)({
  display: 'flex',
  justifyContent: 'space-between',
})

const ContentDelete = styled(Box)`
  position: absolute;
  top: 50%;
  left: 50%;
  transform: translate(-50%, -50%);
  width: 450px;
  background-color: rgb(255, 255, 255);
  box-shadow: rgb(0 0 0 / 20%) 0px 11px 15px -7px,
    rgb(0 0 0 / 14%) 0px 24px 38px 3px, rgb(0 0 0 / 12%) 0px 9px 46px 8px;
  padding: 16px;
  border-radius: 4px;
  outline: none;
`

const BoxConfirm = styled(Box)({
  margin: '20px 0',
})

const FormInline = styled(Box)({
  display: 'flex',
  justifyContent: 'space-between',
  marginBottom: 10,
})

const Label = styled(Typography)({
  fontSize: 14,
  marginTop: 7,
  width: '100%',
})

const ButtonConfirm = styled(Button)({
  height: 36,
  color: '#ffffff',
  marginTop: -1,
  width: 90,
  backgroundColor: '#283237 !important',
  "&:hover": {
    backgroundColor: '#283237',
  }
})

export default ChangePasswordModal
