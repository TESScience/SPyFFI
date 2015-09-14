#!/usr/bin/env bash
screen -S midlatitude -d -m ./midlatitude.py
screen -S sparse -d -m ./sparse.py
screen -S crowded -d -m ./crowded.py
screen -S testpattern -d -m ./testpattern.py
