FROM node:14 as node
COPY . ./
ENV GENERATE_SOURCEMAP false 
RUN yarn install
RUN yarn build

FROM nginx:1.20.1
COPY --from=node ./build /usr/share/nginx/html
CMD nginx -g "daemon off;"