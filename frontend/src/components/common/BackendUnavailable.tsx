import { FC } from "react"

import { Box, Stack, styled, Typography } from "@mui/material"

type Props = {
  // true: background retries are running and the page will auto-recover.
  // false: polling has stopped and the user must reload manually.
  isRetrying: boolean
}

const BackendUnavailable: FC<Props> = ({ isRetrying }) => {
  return (
    <Wrapper>
      <Content>
        <Title>Cannot connect to the OptiNiSt server</Title>
        <Message>
          The server is currently unavailable. Please try again after a while.
        </Message>
        {isRetrying ? (
          <SubMessage>
            This page will automatically recover once the connection is
            restored.
          </SubMessage>
        ) : (
          <SubMessage>Please reload this page to retry.</SubMessage>
        )}
      </Content>
    </Wrapper>
  )
}

const Wrapper = styled(Box)({
  width: "100vw",
  height: "100vh",
  display: "flex",
  alignItems: "center",
  justifyContent: "center",
  padding: 16,
})

const Content = styled(Stack)({
  padding: 32,
  boxShadow: "2px 1px 3px 1px rgba(0,0,0,0.1)",
  borderRadius: 4,
  maxWidth: 480,
  textAlign: "center",
  gap: 12,
})

const Title = styled(Typography)({
  fontSize: 18,
  fontWeight: 600,
})

const Message = styled(Typography)({
  fontSize: 14,
  color: "rgba(0, 0, 0, 0.75)",
})

const SubMessage = styled(Typography)({
  fontSize: 12,
  color: "rgba(0, 0, 0, 0.55)",
})

export default BackendUnavailable
