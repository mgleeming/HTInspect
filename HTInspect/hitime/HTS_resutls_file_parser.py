import os, sys, re

class hitime_hit (object):
    def __init__(self):
        return

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
                hit.rt = rt
                hit.mz = mz
                hit.amp = amp
                print (rt, mz, amp)
                rawData.append(hit)
            else:
                pass

    # ensure that returned data is sorted
    rawData.sort(key = lambda x: x.amp, reverse = True)
    return rawData

