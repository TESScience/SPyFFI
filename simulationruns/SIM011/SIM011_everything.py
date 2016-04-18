from SIM011_base import *

# turn on jitter
inputs['expose']['jitter'] = True
inputs['jitter']['amplifyinterexposurejitter'] = 10.0

# turn on aberration
inputs['camera']['aberrate'] = True
inputs['camera']['warpspaceandtime'] = 0.01

# turn on aberration
inputs['camera']['variablefocus'] = True

# label what kinds of effects are included
inputs['camera']['label'] += '_justaberration'
