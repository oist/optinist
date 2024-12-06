const HOST =
  process.env.REACT_APP_SERVER_HOST || window.location.hostname || "localhost"
const PROTO =
  process.env.REACT_APP_SERVER_PROTO ||
  window.location.protocol.replace(":", "") ||
  "http"
const PORT =
  process.env.REACT_APP_SERVER_PORT ||
  window.location.port ||
  (PROTO == "https" ? 443 : 80)

export const BASE_URL =
  PORT == null ? `${PROTO}://${HOST}` : `${PROTO}://${HOST}:${PORT}`
