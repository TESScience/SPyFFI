from SIM011_base import *

# turn off jitter
inputs['expose']['jitter'] = False
inputs['jitter']['amplifyinterexposurejitter'] = 0.0

# turn off aberration
inputs['camera']['aberrate'] = False

# turn off aberration
inputs['camera']['variablefocus'] = True

# label what kinds of effects are included
inputs['camera']['label'] += '_justfocus'
