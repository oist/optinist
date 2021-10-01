import React from 'react'
// import logo from './logo.svg';

export const ImageViewer = React.memo(() => {
  return <img src={`${process.env.PUBLIC_URL}/sample_image.png`} />
})
