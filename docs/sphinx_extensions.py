import re

from docutils import nodes
from docutils.statemachine import StateMachine
from sphinx.application import Sphinx
from sphinx.directives.other import Include
from sphinx.util.typing import ExtensionMetadata


class IncludeNestml(Include):
    """A directive to extend include for nestml docstrings"""

    def run(self) -> list[nodes.Node]:
        # To properly emit "include-read" events from included RST text,
        # we must patch the ``StateMachine.insert_input()`` method.
        # In the future, docutils will hopefully offer a way for Sphinx
        # to provide the RST parser to use
        # when parsing RST text that comes in via Include directive.
        def _insert_input(include_lines: list[str], source: str) -> None:
            include_lines = [
                re.sub("#[ ]?", "", l_, 1) if l_.startswith("#") else l_ for l_ in include_lines
            ]
            for l_ in include_lines:
                print(l_)

            # Call the parent implementation.
            # Note that this snake does not eat its tail because we patch
            # the *Instance* method and this call is to the *Class* method.
            return StateMachine.insert_input(self.state_machine, include_lines, source)

        self.state_machine.insert_input = _insert_input
        return super().run()


def setup(app: Sphinx) -> ExtensionMetadata:
    app.add_directive("include-nestml", IncludeNestml)

    return {
        "version": "0.1",
        "parallel_read_safe": True,
        "parallel_write_safe": True,
    }
