module.exports = function override(config) {
  okv = {
    fallback: {
      stream: require.resolve('stream-browserify'),
    },
  }
  if (!config.resolve) {
    config.resolve = okv
  } else {
    config.resolve.fallback = { ...config.resolve.fallback, ...okv.fallback }
  }
  return config
}
