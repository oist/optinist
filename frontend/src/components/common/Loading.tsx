import { Box, keyframes, styled } from '@mui/material'

const Loading = () => {
  return (
    <LoaderWrapper>
      <Loader />
    </LoaderWrapper>
  )
}

const LoaderWrapper = styled(Box)(({ theme }) => ({
  top: 0,
  bottom: 0,
  left: 0,
  right: 0,
  position: 'fixed',
  backgroundColor: 'rgba(255,255,255,0.6)',
  zIndex: 100000,
}))

const rotate = keyframes`
  100%   {transform: rotate(360deg)}
`

const prixClipFix = keyframes`
  0%   {clip-path:polygon(50% 50%,0 0,0 0,0 0,0 0,0 0)}
  25%  {clip-path:polygon(50% 50%,0 0,100% 0,100% 0,100% 0,100% 0)}
  50%  {clip-path:polygon(50% 50%,0 0,100% 0,100% 100%,100% 100%,100% 100%)}
  75%  {clip-path:polygon(50% 50%,0 0,100% 0,100% 100%,0 100%,0 100%)}
  100% {clip-path:polygon(50% 50%,0 0,100% 0,100% 100%,0 100%,0 0)}
`

const Loader = styled('span')(({ theme }) => ({
  display: 'block',
  width: 48,
  height: 48,
  borderRadius: '50%',
  position: 'relative',
  zIndex: 100,
  top: 'calc(50% - 24px)',
  left: 'calc(50% - 24px)',
  animation: `${rotate} 1s linear infinite`,
  '&:before': {
    content: "''",
    boxSizing: 'border-box',
    position: 'absolute',
    inset: 0,
    borderRadius: '50%',
    border: '3px solid rgba(0,0,0,0.8)',
    animation: `${prixClipFix} 2s linear infinite`,
  },
}))

export default Loading
