

def get_min_number_of_copies_to_span_read(pattern, read_length=150):
    return int(round(float(read_length) / len(pattern) + 0.499))


def get_gc_content(s):
    res = 0
    for e in s:
        if e == 'G' or e == 'C':
            res += 1
    return float(res) / len(s)
