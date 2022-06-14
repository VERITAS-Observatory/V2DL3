import numpy as np


def decodeConfigMask(mask=15):
    """Decode the telescope config mask to find the telescpes in the array"""
    tels = []
    if mask >= 8:
        tels.append(4)
        mask -= 8
    if mask >= 4:
        tels.append(3)
        mask -= 4
    if mask >= 2:
        tels.append(2)
        mask -= 2
    if mask >= 1:
        tels.append(1)
        mask -= 1
    return sorted(tels)


def produceTelList(mask):
    """Convert the list of telescopes into a string for FITS header"""
    telList = "".join("T" + str(tel) + "," for tel in decodeConfigMask(mask))
    return telList[:-1]


def parseTimeCut(tCutStr):
    """Parse time cuts from time cut text"""
    cut_arr = []
    for cc in tCutStr.split(","):
        start = float(cc.split("/")[0])
        end = float(cc.split("/")[1])
        cut_arr.append((start, end))
    return cut_arr


"""Search a ROOT fCutsFileText for the params named in cut_searches

Returns a dict that will only contain keys for found cuts.
Use the `in` operator on the returned dictionary to see if the cut was found.

Do note that the cut values are saved as strings.

Arguments:
    config_str_ori  -- fCutsFileText
    cut_searches    -- List of the desired parameter names as strings
"""
def getCuts(config_str_ori, cut_searches):
    config_str = str(config_str_ori)
    cuts_found = {}
    for line in config_str.splitlines():
        # Skip comment lines
        if (len(line) == 0) or (line.strip()[0] == "#"):
            continue
        else:
            for cut_search in cut_searches:
                if line.find(cut_search) >= 0:
                    key, cut_str = line.split(" ")
                    if len(cut_str) > 0:
                        cuts_found[cut_search] = cut_str
    return cuts_found


def getTimeCut(config_str_ori):
    """Get time cut from extracted cut config text."""
    config_str = str(config_str_ori)
    for i in config_str.splitlines():
        # Skip comment lines
        if (len(i) == 0) or (i.strip()[0] == "#"):
            continue
        # I
        elif i.find("ES_CutTimes") >= 0:
            key, cut_str = i.split(" ")
            if len(cut_str) == 0:
                return parseTimeCut("0/0")
            return parseTimeCut(cut_str)


def getThetaSquareCut(config_str_ori):
    config_str = str(config_str_ori)
    for i in config_str.splitlines():
        # Skip comment lines
        if (len(i) == 0) or (i.strip()[0] == "#"):
            continue
        # I
        elif i.find("ThetaSquareUpper") >= 0:
            key, cut_str = i.split(" ")
            if len(cut_str) == 0:
                raise Exception("No theta2 cuts present in EA file")
            return float(cut_str)


def isMergable(cut1, cut2):
    """Check if two cuts can be merged."""
    s1, e1 = cut1
    s2, e2 = cut2
    if s1 > s2:
        s1, e1 = cut2
        s2, e2 = cut1
    return s2 <= e1


def mergeTwoTimeCut(cut1, cut2):
    """Merge two cuts."""
    s1, e1 = cut1
    s2, e2 = cut2
    if s1 > s2:
        s1, e1 = cut2
        s2, e2 = cut1
    return s1, max(e1, e2)


def mergeTimeCut(cuts):
    """Go through each cut and merge"""
    cut_sorted = sorted(cuts)
    merged_cut = []
    while len(cut_sorted) > 0:
        test = cut_sorted[0]
        cut_sorted = cut_sorted[1:]
        the_rest = []
        for cc in cut_sorted:
            if isMergable(test, cc):
                test = mergeTwoTimeCut(test, cc)
            else:
                the_rest.append(cc)

        cut_sorted = the_rest
        merged_cut.append(test)

    return list(filter(lambda x: x[0] != x[1], merged_cut))


def getGTArray(startTime_s, endTime_s, cuts):
    """Produce Good Time Interval start and stop time array."""
    if len(cuts) == 0:
        return np.array([startTime_s]), np.array([endTime_s])
    goodTimeStart = [startTime_s] + [cc[1] + startTime_s for cc in cuts]
    goodTimeStop = [cc[0] + startTime_s for cc in cuts] + [endTime_s]

    goodTimeStop[-1] = min(goodTimeStop[-1], endTime_s)
    if goodTimeStart[-1] > goodTimeStop[-1]:
        return goodTimeStart[:-1], goodTimeStop[:-1]
    return np.array(goodTimeStart), np.array(goodTimeStop)
