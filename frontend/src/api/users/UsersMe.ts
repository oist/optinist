import axios from 'utils/axios'
import { UserDTO, UpdateUserDTO, UpdateUserPasswordDTO } from './UsersApiDTO'

export const getMeApi = async (): Promise<UserDTO> => {
  const response = await axios.get('/users/me')
  return response.data
}

export const updateMeApi = async (data: UpdateUserDTO): Promise<UserDTO> => {
  const response = await axios.put('/users/me', data)
  return response.data
}

export const updateMePasswordApi = async (
  data: UpdateUserPasswordDTO,
): Promise<UserDTO> => {
  const response = await axios.put('/users/me/password', data)
  return response.data
}
