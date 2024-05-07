#!/bin/bash

# -u, --update skip files that are newer on the receiver
# -r, --recursive recurse into directories
# -l, --links copy symlinks as symlinks
# -t, --times preserve modification times
# -v, --verbose increase verbosity
# -e, allows you to specify your remote shell, in this case ssh.
# !! --delete is not used here and you should be careful with its usage !!
# $* forwards all arguments that are passed to this script to rsync, e.g.
# --dry-run
rsync -urltv $* -e  ssh . lgold@yael.wagner.lan:rna_folding 
