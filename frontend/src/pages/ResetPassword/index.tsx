import { ChangeEvent, FormEvent, useState } from 'react'
import { useNavigate } from "react-router-dom"
import { Box, Stack, styled, Typography, Link } from '@mui/material'
import { sendResetPasswordMailApi } from 'api/auth/Auth'
import Loading from "components/common/Loading"
import { regexEmail } from 'const/Auth'

const ResetPassword = () => {
    const navigate = useNavigate()
    const [isLoading, setIsLoading] = useState(false)
    const [errors, setErrors] = useState<{ [key: string]: string }>({
        email: '',
    })
    const [values, setValues] = useState<{ email: string }>({
        email: '',
    })

    const onReset = async (event: FormEvent<HTMLFormElement>) => {
        event.preventDefault()

        const errorCheck = validateSubmit()
        if (errors.email || errorCheck) return
        setIsLoading(true)
        try {
            await sendResetPasswordMailApi(values.email)
            setTimeout(()=>{
                alert(` You'll receive a link to reset your password at ${values.email}. Please check your mail!`)
            },300)
        }
        catch {
            setTimeout(()=>{
                alert('Email does not exists!')
            },300)
        }
        finally {
            setIsLoading(false)
        }
    }

    const validateEmail = (value: string): string => {
        if(!value) return 'This field is required'
        if(value.length > 255) return 'The text may not be longer than 255 characters'
        if(!regexEmail.test(value)) return 'The email is invalid'
        return ''
    }

    const validateSubmit = () => {
        let errors = { email: '' }
        errors.email = validateEmail(values.email)
        setErrors(errors)
        return errors.email
    }

    const onChangeValue = (event: ChangeEvent<HTMLInputElement>) => {
        const { name, value } = event.target
        const error = validateEmail(value)
        setValues({ ...values, [name]: value })
        setErrors({ ...errors, [name]: error })
    }

    return (
        <LoginWrapper>
            <LoginContent>
                <Heading>Forgot password?</Heading>
                <Title>No worries, we'll send you reset instructions.</Title>
                <FormSignUp autoComplete="off" onSubmit={onReset}>
                    <Box sx={{ position: 'relative' }}>
                        <LabelField>
                            Email<LableRequired>*</LableRequired>
                        </LabelField>
                        <Input
                            autoComplete="off"
                            error={!!errors.email}
                            name="email"
                            onChange={onChangeValue}
                            value={values.email}
                            placeholder="Enter your email"
                        />
                        <TextError>{errors.email}</TextError>
                    </Box>
                    <Stack
                        flexDirection="row"
                        gap={2}
                        mt={3}
                        alignItems="center"
                        justifyContent="space-between"
                    >
                        <ButtonSignIn onClick={() => navigate('/login')}>Back to SIGN IN</ButtonSignIn>
                        <ButtonLogin type="submit">Reset Password</ButtonLogin>
                    </Stack>
                </FormSignUp>
            </LoginContent>
            {
                isLoading && <Loading />
            }
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

const Heading = styled(Typography)({
    fontSize: 22,
    textAlign: 'center',
    fontWeight: 600,
})

const Title = styled(Typography)({
    textAlign: 'center',
    marginBottom: 24,
    fontSize: 12,
    color: 'rgba(0, 0, 0, 0.65)',
})

const ButtonSignIn = styled(Link)({
    fontSize: 12,
    '&:hover': {
        cursor: 'pointer'
    }
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

export default ResetPassword
