#!/bin/bash


twine upload -u __token__ -p $pypi_token $1
