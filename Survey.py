'''(Maybe), keep track of a TESS survey, with multiple cameras over multiple pointings.'''

from zachopy.Talker import Talker


class Survey(Talker):
    def __init__(self):
        Talker.__init__(self)
        self.offset_from_ecliptic = 6.0									# degrees
        self.number_of_segments = 13
