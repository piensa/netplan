import attr
import datetime
import os
import re
import shlex
import simplejson as json
import stat
from six import string_types
from six.moves.urllib.parse import urlencode as format_query
from subprocess import CalledProcessError, Popen, check_output, PIPE, STDOUT

from .calculator import get_int
from .exceptions import InvisibleRoadsError


@attr.s
class Callback(object):

    id = attr.ib()
    datetime = attr.ib()


def run_command(command, exception_by_error=None):
    if not isinstance(command, string_types):
        command = ' '.join(command)
    command = command.split(';', 1)[0]
    return run_raw_command(command, exception_by_error)


def run_raw_command(command, exception_by_error=None):
    if isinstance(command, string_types):
        command = shlex.split(command)
    try:
        output = check_output(command, stderr=STDOUT)
    except CalledProcessError as e:
        o = e.output
        for error_text, exception in (exception_by_error or {}).items():
            if error_text in o:
                raise exception
        else:
            raise InvisibleRoadsError(o)
    except OSError as e:
        raise InvisibleRoadsError(e.strerror)
    return output.strip()


def schedule_curl_callback(
        minute_count, base_url, value_by_key=None, headers=None, method='GET',
        with_retry=True):
    # Prepare more_lines
    full_url = base_url
    more_lines = []
    if value_by_key:
        if method == 'POST':
            headers['Content-Type'] = 'application/json'
            more_lines.append('-d \'%s\'' % json.dumps(value_by_key))
        else:
            full_url += format_query(value_by_key)
    if with_retry:
        more_lines.extend(['--retry 7'])
    more_lines.extend(['-X %s' % method, full_url])
    # Prepare header_lines
    header_lines = []
    for k, v in (headers or {}).items():
        header_lines.append('-H "%s: %s"' % (k, v))
    # Schedule callback
    shell_parts = ['curl'] + header_lines + more_lines
    shell_text = ' '.join(shell_parts)
    return schedule_shell_callback(minute_count, shell_text)


def schedule_shell_callback(minute_count, shell_text):
    process = Popen([
        'at', 'now + %s minutes' % get_int(minute_count),
    ], stdin=PIPE, stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate((shell_text + '\n').encode('utf-8'))
    try:
        callback_id, when_string = re.search(
            'job (\d+) at (.+)', stderr.decode('utf-8')).groups()
    except AttributeError:
        raise InvisibleRoadsError('cannot schedule callback')
    return Callback(id=int(callback_id), datetime=datetime.datetime.strptime(
        when_string, '%a %b %d %H:%M:%S %Y'))


def cancel_shell_callback(callback_id):
    process = Popen(['atrm', str(callback_id)])
    process.wait()


def format_variables_as_shell_script(d):
    lines = []
    for k, v in d.items():
        if v is None:
            v = ''
        lines.append('%s="%s"' % (k.upper(), v))
    return '\n'.join(lines)


def make_executable(path):
    # https://stackoverflow.com/a/12792002/192092
    x_stat = os.stat(path)
    os.chmod(path, x_stat.st_mode | stat.S_IEXEC)
    return path
