export type UserDTO = {
  uid: string
  email: string
}

export type AddUserDTO = {
  email: string
  password: string
}

export type ListUsersQueryDTO = {
  offset?: number
  limit?: number
}

export type UserListDTO = {
  data: UserDTO[]
  total_page: number
}

export type UpdateUserDTO = {
  email: string
}

export type UpdateUserPasswordDTO = {
  old_password: string
  new_password: string
}
