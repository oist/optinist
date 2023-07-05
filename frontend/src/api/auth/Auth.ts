import axios from 'utils/axios'
import { getRefreshToken } from 'utils/auth/AuthUtils'

export type LoginDTO = {
  email: string
  password: string
}

export type TokenDTO = {
  access_token: string
  token_type: string
  refresh_token: string
  ex_token: string
}
export type AccessTokenDTO = {
  access_token: string
}

export const loginApi = async (data: LoginDTO): Promise<TokenDTO> => {
  const response = await axios.post('/auth/login', data)
  return response.data
}

export const refreshTokenApi = async (): Promise<AccessTokenDTO> => {
  const response = await axios.post('/auth/refresh', {
    refresh_token: getRefreshToken(),
  })
  return response.data
}

export const sendResetPasswordMailApi = async (
  email: string,
): Promise<string> => {
  const response = await axios.post(
    `/auth/send_reset_password_mail?email=${email}`,
  )
  return response.data
}
