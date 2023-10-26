const HOST = process.env.REACT_APP_SERVER_HOST ?? "localhost"
const PORT = process.env.REACT_APP_SERVER_PORT ?? 8000
const PROTO = process.env.REACT_APP_SERVER_PROTO ?? "http"

export const BASE_URL =
  PORT == null ? `${PROTO}://${HOST}` : `${PROTO}://${HOST}:${PORT}`
