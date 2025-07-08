import re

from docutils import nodes
from docutils.statemachine import StateMachine
from sphinx.application import Sphinx
from sphinx.directives.other import Include
from sphinx.util.typing import ExtensionMetadata


class IncludeNestml(Include):
    """A directive to extend include for nestml docstrings"""

    def run(self) -> list[nodes.Node]:
        # To remove the comment symbol at the start of each line in nestml file,
        # we must patch the ``StateMachine.insert_input()`` method.
        def _insert_input(include_lines: list[str], source: str) -> None:
            include_lines = [
                re.sub("#[ ]?", "", l_, 1) if l_.startswith("#") else l_ for l_ in include_lines
            ]
            # Call the parent implementation.
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
