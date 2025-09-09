import sys
import datetime
import traceback
from collections import defaultdict
from typing import DefaultDict

class Logger:
    def __init__(self, level='INFO'):
        self.log_levels = {'INFO': 1, 'DEBUG': 2, 'TRACE': 3}
        self.set_level(level)

    def set_level(self, level):
        if level.upper() not in self.log_levels:
            raise ValueError(f"Unknown logging level: {level}")
        self.current_level = self.log_levels[level.upper()]

    def _log(self, msg, level_name, file=sys.stderr):
        # A private helper method to avoid code repetition
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        formatted_msg = msg.lstrip()
        leading_ws = msg[:len(msg) - len(formatted_msg)]

        print(f"{leading_ws}[{level_name}] [{timestamp}] {formatted_msg}", file=file)

    def info(self, msg):
        if self.current_level >= self.log_levels['INFO']:
            self._log(msg, 'INFO')

    def debug(self, msg):
        if self.current_level >= self.log_levels['DEBUG']:
            self._log(msg, 'DEBUG')

    def trace(self, msg):
        if self.current_level >= self.log_levels['TRACE']:
            self._log(msg, 'TRACE')

    def warning(self, msg):
        # Always print warnings, regardless of level
        self._log(msg, 'WARNING')

    def error(self, msg, exc_info=None):
        """
        Logs an error message, optionally with exception information.
        :param msg: The error message string.
        :param exc_info: An exception object (e.g., from a 'except Exception as e:' block).
        """
        self._log(msg, 'ERROR')
        if exc_info:
            # Format and print the traceback
            tb_str = traceback.format_exc()
            self._log(f"{tb_str}", "ERROR")

_warning_count: DefaultDict[str, int] = defaultdict(int)

def warn_once(logger, msg: str) -> None:
    if _warning_count[msg] == 0 and logger.current_level <= 1:
        logger.warning(msg + " Hiding further warnings of this type. Use logging level DEBUG, or TRACE to see all warnings.")
    else:
        if logger.current_level > 1:
            logger.warning(msg)
    _warning_count[msg] += 1

logger = Logger()