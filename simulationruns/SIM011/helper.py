#!/usr/bin/env python
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='helper for calling SPyFFI with options'
    )
    parser.add_argument('field', type=str, help='what field?')
    parser.add_argument('--jitter',  help='jitter field?', action='store_true')
    parser.add_argument('--aberrate',  help='aberrate positions?', action='store_true')
    parser.add_argument('--focus',  help='variable focus?', action='store_true')
    parser.add_argument('--stamps', help='short cadence postage stamps?', action='store_true')
    args = parser.parse_args()

    from SIM011_base import *

    #inputs['camera']['label'] += '_'
    if 'orion' in args.field.lower():
        # point at orion (crowded)
        inputs['camera']['ra'] = 82.0
        inputs['camera']['dec'] = 1.0
        #inputs['camera']['label'] += 'orion'
    if 'nep' in args.field.lower():
        # point at north ecliptic pole (intermediate)
        inputs['camera']['ra'] = 270.00000
        inputs['camera']['dec'] = 66.56071
        #inputs['camera']['label'] += 'nep'
    if 'ngp' in args.field.lower():
        # point at north galactic pole (sparse)
        inputs['camera']['ra'] = 192.85948
        inputs['camera']['dec'] = 27.12830
        #inputs['camera']['label'] += 'ngp'


    wino = {True:'wi', False:'no'}
    if args.jitter:
        inputs['expose']['jitter'] = True
    else:
        inputs['camera']['jitterkw']['amplifyinterexposurejitter'] = 0.0
    inputs['camera']['label'] += '_{}jitter'.format(wino[args.jitter])

    if args.aberrate:
        inputs['camera']['aberrate'] = True
    inputs['camera']['label'] += '_{}aberrate'.format(wino[args.aberrate])

    if args.aberrate:
        inputs['camera']['variablefocus'] = True
    inputs['camera']['label'] += '_{}focus'.format(wino[args.focus])

    ndays = 27.6
    if args.stamps:
        inputs['camera']['stamps'] = 4000
        inputs['observation']['cadencestodo'] = {120:int(ndays*48*15)}
    else:
        inputs['camera']['stamps'] = None
        inputs['observation']['cadencestodo'] = {1800:int(ndays*48)}

    print inputs['camera']
    print inputs['observation']

    o = Observation(inputs)
    o.create()
