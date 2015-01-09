class Talker(object):
    '''Objects the inherit from Talker have "mute" and "pithy" attributes, a report('uh-oh!') method that prints when unmuted, and a speak('yo!') method that prints only when unmuted and unpithy.'''
    def __init__(self, mute=False, pithy=False):
        self.mute = False
        self.pithy = False

    def speak(self, string, level=0):
        '''If verbose=True and terse=False, this will print to terminal. Otherwise, it won't.'''
        if self.pithy == False:
            self.report(string, level)

    def report(self, string, level=0):
        '''If verbose=True, this will print to terminal. Otherwise, it won't.'''
        if self.mute == False:
            self.prefix = '{spacing}[{name}] '.format(name = self.__class__.__name__.lower(), spacing = ' '*level)
            equalspaces = ' '*len(self.prefix)
            print self.prefix + string.replace('\n', '\n' + equalspaces)
