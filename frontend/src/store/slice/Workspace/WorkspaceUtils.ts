import { UserDTO } from "api/users/UsersApiDTO"

export const isMine = (user?: UserDTO, idUserWorkSpace?: number) => {
  return !!(user && idUserWorkSpace && user.id === idUserWorkSpace)
}
