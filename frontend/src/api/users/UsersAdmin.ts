import axios from 'utils/axios'
import {
  AddUserDTO,
  UserDTO,
  ListUsersQueryDTO,
  UserListDTO,
  UpdateUserDTO,
} from './UsersApiDTO'

export const createUserApi = async (data: AddUserDTO): Promise<UserDTO> => {
  const response = await axios.post('/admin/users', data)
  return response.data
}

export const getUserApi = async (uid: string): Promise<UserDTO> => {
  const response = await axios.get(`/admin/users/${uid}`)
  return response.data
}

export const listUsersApi = async (
  data: ListUsersQueryDTO,
): Promise<UserListDTO> => {
  const response = await axios.get('/admin/users', { params: data })
  return response.data
}

export const updateUserApi = async (
  uid: string,
  data: UpdateUserDTO,
): Promise<UserDTO> => {
  const response = await axios.put(`/admin/users/${uid}`, data)
  return response.data
}

export const deleteUserApi = async (uid: string): Promise<string> => {
  const response = await axios.delete(`/admin/users/${uid}`)
  return response.data
}
