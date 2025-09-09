# Code taken from WhatsHap (https://github.com/whatshap/whatshap)

import sys
import pkgutil
import importlib

import giggles.cli as cli_package
from giggles.logger import logger
from . import __version__
from .args import HelpfulArgumentParser


def ensure_pysam_version():
    from pysam import __version__ as pysam_version
    from distutils.version import LooseVersion

    if LooseVersion(pysam_version) < LooseVersion("0.8.1"):
        sys.exit("Giggles requires pysam >= 0.8.1")


def main(argv=sys.argv[1:]):
    ensure_pysam_version()
    parser = HelpfulArgumentParser(description=__doc__, prog="giggles")
    parser.add_argument("--version", action="version", version="%(prog)s " + __version__)
    parser.add_argument("--logging-level", default="INFO", choices=["INFO", "DEBUG", "TRACE"], help="Set the logging level. (Level: INFO < DEBUG < TRACE)")
    subparsers = parser.add_subparsers()

    # Import each module that implements a subcommand and add a subparser for it.
    # Each subcommand is implemented as a module in the cli subpackage.
    # It needs to implement an add_arguments() and a main() function.
    modules = pkgutil.iter_modules(cli_package.__path__)
    for _, module_name, _ in modules:
        module = importlib.import_module("." + module_name, cli_package.__name__)
        subparser = subparsers.add_parser(
            module_name,
            help=module.__doc__.strip().split("\n", maxsplit=1)[0],
            description=module.__doc__,
        )
        subparser.set_defaults(module=module, subparser=subparser)
        module.add_arguments(subparser)

    args = parser.parse_args(argv)
    logger.set_level(args.logging_level)

    if not hasattr(args, "module"):
        parser.error("Please provide the name of a subcommand to run")
    else:
        module = args.module
        if hasattr(args.module, "validate"):
            subparser = args.subparser
            args.module.validate(args, subparser)
        del args.subparser
        del args.module
        del args.logging_level
        try:
            module.main(args)
        except Exception as e:
            logger.error(f"giggles error: {str(e)}", exc_info=e)
            sys.exit(1)


if __name__ == "__main__":
    main()
