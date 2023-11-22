import axios from "utils/axios"

export const getModeStandaloneApi = async (): Promise<boolean> => {
  const res = await axios.get("/is_standalone")
  if (res.data === undefined) {
    return getModeStandaloneApi()
  }
  return !!res.data
}
