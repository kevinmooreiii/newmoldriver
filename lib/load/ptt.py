""" Library of patterns to simplify the parsing of input files
"""

import autoparse.find as apf
import autoparse.pattern as app


# Patterns that will be helpful
KEYWORD_KEYVALUE_PATTERN = (
    app.capturing(app.one_or_more(app.NONSPACE)) +
    app.zero_or_more(app.SPACE) +
    '=' +
    app.zero_or_more(app.SPACE) +
    app.capturing(app.one_or_more(app.NONSPACE)))


def read_inp_str(filepath, filename):
    """ read the run parameters from a file
    """
    input_file = os.path.join(filepath, filename)
    with open(input_file, 'r') as inp_file:
        inp_lines = inp_file.readlines()
    inp_str = ''.join(
        (line for line in inp_lines
         if '#' not in line and line.strip() != ''))
    return inp_str


# Reading various section strings from the files #
def paren_section(string):
    """ Read the string that has the global model information
    """
    return (string + app.SPACE + app.escape('=') + app.SPACE +
            app.escape('(') +
            app.capturing(app.one_or_more(app.WILDCARD, greedy=False)) +
            app.escape(')'))


def end_section(string):
    """ Read the string that has the global model information
    """
    return (string + app.one_or_more(app.SPACE) +
            app.capturing(app.one_or_more(app.WILDCARD, greedy=False)) +
            'end')


def end_section_wname(string, name):
    """ Read the string that has the global model information
    """
    return (string + app.one_or_more(app.SPACE) + name +
            app.capturing(app.one_or_more(app.WILDCARD, greedy=False)) +
            'end')


def end_section_wname2(string):
    """ Read the string that has the global model information
    """
    return (string + app.SPACES +
            app.capturing(app.one_or_more(app.NONSPACE)) +
            app.capturing(app.one_or_more(app.WILDCARD, greedy=False)) +
            'end')


def keyword_pattern(string):
    """ Generates the key pattern string
    """
    return (string +
            app.zero_or_more(app.SPACE) +
            '=' +
            app.zero_or_more(app.SPACE) +
            app.capturing(app.one_or_more(app.NONSPACE)))


def parse_idx_inp(idx_str):
    """ parse idx string
    """
    idx_str = idx_str.strip()
    if idx_str.isdigit():
        idxs = [int(idx_str), int(idx_str)]
    if '-' in idx_str:
        [idx_begin, idx_end] = idx_str.split('-')
        idxs = list(range(int(idx_begin), int(idx_end)+1))
    return idxs


def remove_comment_lines(section_str):
    """ Remove lines of comments from strings of a section
    """
    section_lines = section_str.splitlines()
    cleaned_str = ''.join([line for line in section_lines
                           if '#' not in line])
    return cleaned_str


# Build a keyword dictionary
def build_keyword_dct(section_str):
    """ Take a section with keywords defined and build
        a dictionary for the keywords
    """
    keyword_dct = {}
    for line in section_str.splitlines():
        # Put a cleaner somehwere to get rid of blank lines
        if line.strip() != '':
            [word, val] = apf.first_capture(KEYWORD_KEYVALUE_PATTERN, line)
            keyword_dct[word] = val
    return keyword_dct


# Helper functions
def remove_empty_lines(string):
    """ Remove any empty lines from a string
    """
    return '\n'.join([line for line in string.splitlines()
                      if line.strip()])


# def get_key_str(inp_str, lvl, key):
#     """ Find the value for a given key
#     """
#     key_str = apf.first_capture(
#         keyword_pattern(key), get_lvl_str(inp_str, lvl))
#     assert key_str is not None
#     return key_str
#
#
# def get_key_int(inp_str, lvl, key):
#     """ Find the value for a given key; make integer
#     """
#     key_str = apf.first_capture(
#         keyword_pattern(key), get_lvl_str(inp_str, lvl))
#     assert key_str is not None
#     return int(key_str)
#
#
# Parsing lines within the various section strings
# def channel_block(inp_str):
#     """ get the channels
#     """
#     # Various Patterns for entering the Channel indices
#     cidx_ptt = one_of_these([
#         INTEGER,
#         INTEGER + '-' + INTEGER,
#         series(INTEGER, ',')
#     ])
#     model_pattern = ('channels' + app.one_or_more(app.SPACE) +
#                      app.capturing(cidx_ptt) +
#                      app.zero_or_more(app.SPACE) + NEWLINE +
#                 app.capturing(app.one_or_more(app.WILDCARD, greedy=False)) +
#                      'end')
#     channel_str = apf.first_capture(model_pattern, inp_str)
#
#     # Parse the channel indices into a list of integers for all channels
#     cidx_str = channel_str[0]
#     if ',' in cidx_str:
#         cidxs = cidx_str.strip().split(',')
#     if '-' in cidx_str:
#         [idx_begin, idx_end] = cidx_str.strip().split('-')
#         cidxs = list(range(int(idx_begin), int(idx_end)+1))
#
#     # Get the block of text containing the keyword
#     channel_keys = channel_str[1]
#
#     return cidxs, channel_keys
