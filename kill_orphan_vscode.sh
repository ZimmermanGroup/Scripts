#!/bin/bash

ps -ef | grep $USER | grep "vscode-server" > kill_output
cat kill_output | awk ' { print $2 } ' > kill_list
for i in $(cat kill_list )
do
  echo "kill -9 $i"
  kill -9 $i
done

rm -f kill_output
rm -f kill_list
