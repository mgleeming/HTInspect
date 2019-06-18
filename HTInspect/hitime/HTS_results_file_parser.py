import os, sys, re

class hitime_hit (object):
    def __init__(self):
        return

class raw_data (object):
    def __init__(self):
        self.mz = []
        self.rt = []
        self.amp = []

def parse_line(line):
    numbers = re.findall(r"[-+]?\d*\.\d+|\d+", line)
    if len(numbers) > 0:
        return numbers
    else:
        return None


def reader(inFile):

    rawData = []
    maxScore = 0
    with open(inFile,'r') as data:
        for line in data:
            datum = parse_line(line)
            if datum:
                datum = [float(x) for x in datum]
                try:
                    rt, mz, amp = datum
                except:
                    print ('Error reading line:')
                    print (line)
                    continue

                hit = hitime_hit()
                hit.rt = rt/60
                hit.mz = mz
                hit.amp = amp
                rawData.append(hit)
            else:
                pass

    # ensure that returned data is sorted
    rawData.sort(key = lambda x: x.amp, reverse = True)
    return rawData

def processed_results_reader(if1):

    # flow controls
    read_raw_data = 0
    read_hit_data = 0
    read_headers = 0
    count_hits = 1

    # init raw data structure
    rawData = raw_data()
    HT_hits = []
    headers = {}

    hit = None

    # read file
    with open(if1,'r') as data:
        for i, line in enumerate(data):
            l = line.strip()
            if 'BEGIN_POSTPROCESSING_PARAMETERS' in l: read_headers = 1; continue
            if 'END_POSTPROCESSING_PARAMETERS' in l: read_headers = 0; continue
            if 'BEGIN_RAW_HITIME_RESULTS' in l: read_raw_data = 1; continue
            if 'END_RAW_HITIME_RESULTS' in l: read_raw_data = 0; continue
            if 'BEGIN_Results_for_hit' in l: read_hit_data = 1; hit = hitime_hit(); continue
            if 'END_Results_for_hit' in l: read_hit_data = 0; continue

            if read_headers == 1:
                param, value = l.replace(' ','').split(':::')
                try: value = float(value)
                except: pass
                headers[param.strip('#').strip(' ')] = value

            if read_raw_data == 1:
                datum = parse_line(l)
                if datum is not None and '##' not in line:
                    rt, mz, amp = datum
                    rawData.rt.append(rt)
                    rawData.mz.append(mz)
                    rawData.amp.append(amp)

            if read_hit_data == 1:
                    if '<RT>' in line: hit.rt = float(parse_line(l)[0])
                    if '<MZ>' in line: hit.mz = float(parse_line(l)[0])
                    if '<SCORE>' in line: hit.score = float(parse_line(l)[0])
                    if '<EIC_RT>' in line: hit.EIC_RT = [float(x) for x in parse_line(l)]
                    if '<EIC_int_light>' in line: hit.EIC_int_light = [float(x) for x in parse_line(l)]
                    if '<EIC_int_heavy>' in line: hit.EIC_int_heavy = [float(x) for x in parse_line(l)]
                    if '<MS_mz>' in line: hit.MS_mz = [float(x) for x in parse_line(l)]
                    if '<MS_int>' in line:
                        hit.MS_int = [float(x) for x in parse_line(l)]
                        hit.index = count_hits
                        count_hits += 1
                        HT_hits.append(hit)

    return rawData, HT_hits, headers

