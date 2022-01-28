import os


class RunlistValidationError(Exception):
    pass


class RunlistParsingError(Exception):
    pass


def parseSectionTag(l):
    isEnd = False
    if l[-1] != ']':
        raise RunlistParsingError('Tag not ended with ]')
    l_content = l.strip('[').strip(']')
    if l_content[0] == '/':
        isEnd = True
        l_content.strip('/')
    contents = l_content.split()
    if (len(contents) < 3) or (contents[1] != 'ID:'):
        raise RunlistParsingError('Wrong tag format.')
    key, id_str, gid = contents
    try:
        gid = int(gid)
    except ValueError:
        raise RunlistParsingError(f'Group ID {gid} is not an integer')
    return isEnd, key, gid


def validateRunlist(r_dict):
    must_have = ['EA', 'RUNLIST']
    # Check if EA and RUNLIST are present
    for k in must_have:
        if k not in r_dict.keys():
            raise RunlistValidationError('Missing {}'.format(k))
    # Cehck if there's only one EA per group
    for gid, item in r_dict['EA'].items():
        if len(item) != 1:
            exception_str = 'Only one EA file can be used for each group.\n'
            exception_str += '[EA ID: {:d}] have {:d} files.'.format(gid, len(item))
            raise RunlistValidationError(exception_str)
    # Check if all files exists
    for gid, item in r_dict['EA'].items():
        if not os.path.exists(item[0]):
            raise RunlistValidationError('{} does not exist.'.format(item[0]))

    for gid, item in r_dict['RUNLIST'].items():
        for f in item:
            if not os.path.exists(f):
                raise RunlistValidationError('{} does not exist.'.format(f))

    # Check if all groups have corresponding EA
    rl_minus_ea = list(set(r_dict['RUNLIST'].keys()) - set(r_dict['EA'].keys()))
    if rl_minus_ea:
        raise RunlistValidationError('ID: {} have no matching EA.'.format(rl_minus_ea))


def parseRunlistStrs(lines):
    not_used_str = []
    tag_region_content = []
    lines_sans_whitespace = map(lambda x: x.strip(), lines)
    in_tag_region = False

    current_key_id = None

    encounter_first_tag = False
    parse_dict = {}
    for l in lines_sans_whitespace:
        if (len(l) == 0) or l[0] == '#':
            # Skip empty lines
            continue
        elif l[0] == '[':
            encounter_first_tag = True
            isEnd, key, gid = parseSectionTag(l)
            if (not in_tag_region) and (isEnd):
                raise RunlistParsingError('No start tag found for [{} ID: {:d}]'.format(key, gid))
            elif in_tag_region and (not isEnd):
                raise RunlistParsingError('No ending tag found for previous block before starting [{} ID: {:d}].'.format(key, gid))
            elif not isEnd:
                current_key_id = (key, gid)
                try:
                    parse_dict[key][gid] = []
                except KeyError:
                    parse_dict[key] = {gid: []}

                in_tag_region = True
            else:
                in_tag_region = False
                current_key_id = None
        else:
            # Fill group 0 files
            if not (in_tag_region or encounter_first_tag):
                not_used_str.append(l)
            if in_tag_region:
                key, gid = current_key_id
                parse_dict[key][gid].append(l)
    if in_tag_region:
        raise RunlistParsingError('Missing closing tags.')
    if not_used_str:
        try:
            parse_dict['RUNLIST'][0] = not_used_str
        except KeyError:
            parse_dict['RUNLIST'] = {0: not_used_str}
    return parse_dict
