#!/usr/bin/env bash
screen -S midlatitude -d -m ./midlatitude.py
screen -S sparse -d -m ./sparse.py
screen -S crowded -d -m ./crowded.py
screen -S littletestpattern -d -m ./littletestpattern.py
#screen -S bigtestpattern -d -m ./bigtestpattern.py
