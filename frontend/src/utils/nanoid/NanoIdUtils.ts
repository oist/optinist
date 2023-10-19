import { NANO_ID_LENGTH } from "const/flowchart"
import { customAlphabet } from "nanoid/non-secure"

const nanoid = customAlphabet(
  "0123456789abcdefghijklmnopqrstuvwxyz",
  NANO_ID_LENGTH,
)

export const getNanoId = () => {
  return nanoid()
}
