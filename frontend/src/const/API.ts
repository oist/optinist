const HOST = process.env.REACT_APP_SERVER_HOST
const PORT = process.env.REACT_APP_SERVER_PORT

export const BASE_URL = `http://${HOST}:${PORT}`
export const WS_BASE_URL = `ws://${HOST}:${PORT}`
