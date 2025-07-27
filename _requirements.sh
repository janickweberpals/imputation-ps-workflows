#!/bin/sh

R -e "install.packages('pak')"

R -e "pak::pak('V8')"
