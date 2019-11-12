#!/bin/bash
# script to read submission command from a parameter file

if [ "$1" = "-h" ]; then
# begin help message
echo "
readSubmissionCommand.sh [file with commands]

optional parameters:

    [file with commands]    file containing submission command for scripts
                            (default: submissionCommands.dat)

--------------------------------------------------------------------------------
"
#end help message
exit
fi

# Set submission commands file to arg or default
[[ "$1" ]] && CMDFILE=$1 || CMDFILE="$EVNDISPSYS/scripts/VTS/submissionCommands.dat"

# Check to make sure submission commands file exists
if [[ ! -f "$CMDFILE" ]]; then
    echo "ERROR! Submission command list $CMDFILE not found, exiting..."
    exit 1
fi

# Use submission command with * in front; if multiple commands have an *,
# the command that is farthest down in the file will be used
while read STAR LINE; do
    if [[ $STAR = "*" ]]; then
        CMD=$LINE
    fi
done < $CMDFILE

if [[ -z "$CMD" ]]; then
    echo "ERROR! No submission command is selected in $CMDFILE."
    exit 1
else 
    echo $CMD
fi

exit
