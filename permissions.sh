#!/bin/bash
# set file permisions recursively

setfacl -R -m "g:psm:rwx" $1
setfacl -R -m "g:plb:rwx" $1
setfacl -R -m "u:wangha73:rwx" $1
