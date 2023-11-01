from six.moves.urllib_parse import urlparse as parse_url


PRONOUNCEABLE_ALPHABET = '23456789abcdefghijkmnpqrstuvwxyz'


def format_url(origin='127.0.0.1', port=None, scheme=None):
    if not origin:
        origin = '127.0.0.1'
    if not origin.startswith('http'):
        origin = 'http://' + origin
    x = parse_url(origin)
    port = port or x.port
    scheme = scheme or x.scheme
    if scheme == 'https' and (not port or port == 443):
        url_template = '{scheme}://{hostname}'
        port = 443
    elif scheme == 'http' and (not port or port == 80):
        url_template = '{scheme}://{hostname}'
        port = 80
    else:
        url_template = '{scheme}://{hostname}:{port}'
    return url_template.format(scheme=scheme, hostname=x.hostname, port=port)


def normalize_url(x):
    x = x.strip().lower()
    x = x.replace(' ', '-')
    x = x.replace('_', '-')
    return x


def encode_number(non_negative_integer, alphabet=PRONOUNCEABLE_ALPHABET):
    # http://stackoverflow.com/a/1119769/192092
    if non_negative_integer < 0:
        raise ValueError
    if non_negative_integer == 0:
        return alphabet[0]
    characters = []
    base = len(alphabet)
    while non_negative_integer:
        remainder = non_negative_integer % base
        non_negative_integer = non_negative_integer // base
        characters.append(alphabet[remainder])
    characters.reverse()
    return ''.join(characters)


def decode_number(string, alphabet=PRONOUNCEABLE_ALPHABET):
    # http://stackoverflow.com/a/1119769/192092
    base = len(alphabet)
    string_length = len(string)
    number = 0
    try:
        for character_number, character in enumerate(string, 1):
            power = string_length - character_number
            number += alphabet.index(character) * (base ** power)
    except ValueError:
        number = -1
    return number
