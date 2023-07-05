import { Box, Stack, styled, Typography } from '@mui/material'
import { useDispatch } from 'react-redux'
import { getMe, login } from 'store/slice/User/UserActions'
import { AppDispatch } from 'store/store'
import { ChangeEvent, FormEvent, useState } from 'react'
import { Link, useNavigate } from 'react-router-dom'
import Loading from 'components/common/Loading'

const Login = () => {
  const navigate = useNavigate()
  const dispatch: AppDispatch = useDispatch()

  const [isLoading, setIsLoading] = useState(false)
  const [errors, setErrors] = useState<{ [key: string]: string }>({
    email: '',
    password: '',
  })
  const [values, setValues] = useState<{ email: string; password: string }>({
    email: '',
    password: '',
  })

  const onSubmit = async (event: FormEvent<HTMLFormElement>) => {
    event.preventDefault()
    const errorCheck = validateSubmit()
    if (errorCheck) return
    setIsLoading(true)

    try {
      dispatch(login(values)).unwrap().then((_) => {
        dispatch(getMe()).unwrap().then((_) => {
          navigate('/')
        })
      })
    } catch (e) {
      setErrors({ email: 'Email or password is wrong', password: '' })
    } finally {
      setIsLoading(false)
    }
  }


  const validateSubmit = () => {
    let errors = { email: '', password: '' }
    if (!values.email) {
      errors.email = 'This field is required'
    }
    if (!values.password) {
      errors.password = 'This field is required'
    }
    setErrors(errors)
    return errors.password || errors.email
  }

  const onChangeValue = (event: ChangeEvent<HTMLInputElement>) => {
    const { name, value } = event.target
    setValues({ ...values, [name]: value })
    setErrors({ ...errors, [name]: !value ? 'This field is required' : '' })
  }

  return (
    <LoginWrapper>
      <LoginContent>
        <Title data-testid="title">Sign in to your account</Title>
        <FormSignUp autoComplete="off" onSubmit={onSubmit}>
          <Box sx={{ position: 'relative' }}>
            <LabelField>
              Email<LableRequired>*</LableRequired>
            </LabelField>
            <Input
              data-testid="email"
              autoComplete="off"
              error={!!errors.email}
              name="email"
              onChange={onChangeValue}
              value={values.email}
              placeholder="Enter your email"
            />
            <TextError data-testid="error-email">{errors.email}</TextError>
          </Box>
          <Box sx={{ position: 'relative' }}>
            <LabelField>
              Password<LableRequired>*</LableRequired>
            </LabelField>
            <Input
              data-testid="password"
              autoComplete="off"
              error={!!errors.password}
              onChange={onChangeValue}
              name="password"
              type="password"
              value={values.password}
              placeholder="Enter your password"
            />
            <TextError data-testid="error-password">
              {errors.password}
            </TextError>
          </Box>
          <Description>
            Forgot your password?
            <LinkWrappper to="/reset-password">Reset password</LinkWrappper>
          </Description>
          <Stack
            flexDirection="row"
            gap={2}
            mt={3}
            alignItems="center"
            justifyContent="flex-end"
          >
            <ButtonLogin data-testid="button-submit" type="submit">
              SIGN IN
            </ButtonLogin>
          </Stack>
        </FormSignUp>
      </LoginContent>
      {isLoading && <Loading />}
    </LoginWrapper>
  )
}

const LoginWrapper = styled(Box)({
  width: '100%',
  height: '100%',
  display: 'flex',
  alignItems: 'center',
  justifyContent: 'center',
})

const LoginContent = styled(Box)({
  padding: 30,
  boxShadow: '2px 1px 3px 1px rgba(0,0,0,0.1)',
  borderRadius: 4,
})

const Title = styled(Typography)({
  fontSize: 15,
  fontWeight: 600,
  marginBottom: 24,
})

const FormSignUp = styled('form')({})

const LabelField = styled(Typography)({
  fontSize: 14,
})

const LableRequired = styled('span')({
  color: 'red',
  fontSize: 14,
  marginLeft: 2,
})

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
    marginBottom: 22,
    transition: 'all 0.3s',
    outline: 'none',
    ':focus, :hover': {
      borderColor: '#1677ff',
    },
  }
})

const Description = styled(Typography)(({ theme }) => ({
  fontSize: 12,
  color: 'rgba(0, 0, 0, 0.65)',
  marginTop: theme.spacing(1),
}))

const LinkWrappper = styled(Link)({
  marginLeft: 6,
  color: '#1892d1',
})

const ButtonLogin = styled('button')({
  backgroundColor: '#283237',
  color: '#ffffff',
  borderRadius: 4,
  border: 'none',
  outline: 'none',
  padding: '10px 20px',
  cursor: 'pointer',
})

const TextError = styled(Typography)({
  fontSize: 12,
  color: 'red',
  position: 'absolute',
  bottom: 4,
})

export default Login
