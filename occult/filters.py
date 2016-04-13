from collections import namedtuple
Filter = namedtuple('Filter', 'pivot fwhm')

sdss_filters = {'u': Filter(3540, 570), 'g': Filter(4770, 1370),
                'r': Filter(6230, 1370), 'i': Filter(7630, 1530),
                'z': Filter(9130, 950)}