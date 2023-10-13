import { customAlphabet } from 'nanoid/non-secure'

const nanoid = customAlphabet('0123456789abcdefghijklmnopqrstuvwxyz', 8)

export const getNanoId = () => {
  return nanoid()
}
