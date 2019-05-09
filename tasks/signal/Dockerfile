FROM alpine:3.8

RUN mkdir -p /app
COPY src/ /app/src
COPY Makefile /app

RUN apk add --no-cache build-base jsoncpp-dev zlib-dev cmake wget git lapack-dev curl-dev libgomp && \
    git clone https://github.com/pezmaster31/bamtools && cd bamtools && mkdir build && \
    cd build && cmake .. && make && make install && cd ../.. && rm -r bamtools && \
    git clone https://gitlab.com/conradsnicta/armadillo-code && cd armadillo-code && \
    cmake . && make && make install && cd .. && rm -r armadillo-code
RUN cd /app && make && mv bin/* /bin && apk del build-base cmake wget git
