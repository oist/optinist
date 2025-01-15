/**
 * Collect system version information from package json and export to temporary file
 */

const fs = require("fs")
const packageJson = require("./package.json")

const versionsJson = {
  version: packageJson.version,
  nodeVersion: process.version,
  reactVersion: packageJson.dependencies.react,
}

fs.writeFileSync("./src/.versions.json", JSON.stringify(versionsJson))
