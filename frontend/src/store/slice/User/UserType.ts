import { UserListDTO, UserDTO } from 'api/users/UsersApiDTO'

export const USER_SLICE_NAME = 'user'

export type User = {
  currentUser?: UserDTO
  listUserSearch?: UserDTO[]
  loading: boolean
  listUser?: UserListDTO
}
