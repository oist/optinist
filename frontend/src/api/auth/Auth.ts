import axios from 'utils/axios'
import { getRefreshToken } from 'utils/auth/AuthUtils'

export const loginApi = async (data: { email: string; password: string }) => {
  const response = await axios.post('/auth/login', data)
  return response.data
}

export const refreshTokenApi = async () => {
  const response = await axios.post('/auth/refresh', {
    refresh_token: getRefreshToken(),
  })
  return response.data
}

export const sendResetPasswordMailApi = async (email: string) => {
  const response = await axios.post(
    `/auth/send_reset_password_mail?email=${email}`,
  )
  return response.data
}
