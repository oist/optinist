FROM node:16 as node
COPY . ./
RUN npm install --legacy-peer-deps && npm run build

FROM nginx:1.20.1
COPY --from=node ./build /usr/share/nginx/html
CMD nginx -g "daemon off;"