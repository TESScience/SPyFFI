
for field in ['orion ', 'nep ', 'ngp ']:
    for stamp in ['--stamps ', '']:
        for options in ['--jitter --aberrate --focus',
                        '--jitter',
                        '--aberrate',
                        '--focus',
                        '']:
            command = field + stamp + options
            sessionname = command.replace('--','').replace(' ', '_')
            print 'screen -dmS {} python helper.py {}'.format(sessionname, command)
