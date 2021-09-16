#!/bin/bash
# This is a comment!
echo Running interact script
interact -c $1 -p shared -t $2:00:00
exec bash
