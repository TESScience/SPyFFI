from SIM011_base import *

# turn on jitter
inputs['expose']['jitter'] = True
inputs['jitter']['amplifyinterexposurejitter'] = 10.0

# turn off aberration
inputs['camera']['aberrate'] = False

# turn off aberration
inputs['camera']['variablefocus'] = False

# label what kinds of effects are included
inputs['camera']['label'] += '_justjitter'
