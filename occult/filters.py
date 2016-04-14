from collections import namedtuple
Filter = namedtuple('Filter', 'name pivot fwhm')

sdss_filters = {'u': Filter('u', 3540, 570), 'g': Filter('g', 4770, 1370),
                'r': Filter('r', 6230, 1370), 'i': Filter('i', 7630, 1530),
                'z': Filter('z', 9130, 950)}
                
ha = Filter('ha', 6565, 50)