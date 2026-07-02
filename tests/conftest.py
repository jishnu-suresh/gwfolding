"""Test configuration: inject the fake ``Fr`` module before gwfolding imports it.

gwfolding imports ``Fr`` lazily (inside functions), so registering the stub in
``sys.modules`` here makes every ``import Fr`` / ``from Fr import ...`` resolve to
our synthetic implementation during the test session.
"""

import sys
import os

sys.path.insert(0, os.path.dirname(__file__))

import frstub  # noqa: E402

sys.modules["Fr"] = frstub
