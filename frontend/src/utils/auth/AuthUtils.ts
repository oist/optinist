export const saveToken = (access_token: string) => {
  localStorage.setItem('access_token', access_token)
}

export const getToken = () => {
  return localStorage.getItem('access_token')
}

export const saveRefreshToken = (refresh_token: string) => {
  localStorage.setItem('refresh_token', refresh_token)
}

export const getRefreshToken = () => {
  return localStorage.getItem('refresh_token')
}

export const removeRefreshToken = () => {
  return localStorage.removeItem('refresh_token')
}

export const logout = () => {
  removeRefreshToken()
  removeToken()
  removeExToken()
  window.location.href = '/login'
}

export const removeToken = () => {
  return localStorage.removeItem('access_token')
}

export const saveExToken = (ExToken: string) => {
  localStorage.setItem('ExToken', ExToken)
}

export const getExToken = () => {
  return localStorage.getItem('ExToken')
}

export const removeExToken = () => {
  return localStorage.removeItem('ExToken')
}
