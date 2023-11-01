class InvisibleRoadsError(Exception):
    pass


class BadArchive(IOError, InvisibleRoadsError):
    pass


class BadFormat(IOError, InvisibleRoadsError):
    pass


class BadPath(IOError, InvisibleRoadsError):
    pass
