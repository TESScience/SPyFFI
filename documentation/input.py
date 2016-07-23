documented_keywords = {
    'catalog': {'name': 'Name of catalog to use',
                'skykw': '???',
                'lckw': '???'},
    'observation': {'cadencestodo': '???',
                    'collate': '???'},
    'camera': {'cadence': 'Do not use, use observation.cadencestodo',
               'ra': 'Right Ascension of field center',
               'dec': 'Declination of field center',
               'testpattern': 'Broken, do not use',
               'subarray': 'None for the official 2048**2 pixels, else an integer n for n**2 pixels',
               'label': 'Part of the name of the directory that SPyFFI will create for the output',
               'cameranumber': 'Unimplemented, do not use',
               'warpspaceandtime': 'If a number, use it to scale the speed of light in the aberration calculation. '
                                   'Otherwise, use False',
               'abberate': 'If True, apply abberation of light to star positions',
               'counterstep': 'Normally 1. Scales the time between exposures',
               'positionangle': 'Unimplemented, do not use',

               },
    'expose': {},
    'jitter': {}
}
