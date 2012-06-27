from model import _Call

cdef _map(func, iterable, bad='.'):
    '''``map``, but make bad values None.'''
    return [func(x) if x != bad else None
            for x in iterable]

cdef char *INTEGER = 'Integer'
cdef char *FLOAT = 'Float'
cdef char *NUMERIC = 'Numeric'

def parse_samples(
        list names, list samples, list samp_fmt,
        list samp_fmt_types, list samp_fmt_nums, site):

    cdef char *name, *fmt, *entry_type, *sample
    cdef int i, j
    cdef list samp_data = []
    cdef dict sampdict
    cdef list sampvals
    n_samples = len(samples)
    n_formats = len(samp_fmt)

    for i in range(n_samples):
        name = names[i]
        sample = samples[i]

        # parse the data for this sample
        sampdict = dict([(x, None) for x in samp_fmt])

        sampvals = sample.split(':')

        for j in range(n_formats):
            if j >= len(sampvals):
                break
            fmt = samp_fmt[j]
            vals = sampvals[j]
            entry_type = samp_fmt_types[j]
            # TODO: entry_num is None for unbounded lists
            entry_num = samp_fmt_nums[j]

            # short circuit the most common
            if vals == '.' or vals == './.':
                sampdict[fmt] = None
                continue

            # we don't need to split single entries
            if entry_num == 1 or ',' not in vals:

                if entry_type == INTEGER:
                    sampdict[fmt] = int(vals)
                elif entry_type == FLOAT or entry_type == NUMERIC:
                    sampdict[fmt] = float(vals)
                else:
                    sampdict[fmt] = vals

                if entry_num != 1:
                    sampdict[fmt] = (sampdict[fmt])

                continue

            vals = vals.split(',')

            if entry_type == INTEGER:
                sampdict[fmt] = _map(int, vals)
            elif entry_type == FLOAT or entry_type == NUMERIC:
                sampdict[fmt] = _map(float, vals)
            else:
                sampdict[fmt] = vals

        # create a call object
        call = _Call(site, name, sampdict)
        samp_data.append(call)

    return samp_data
