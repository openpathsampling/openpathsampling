import sys
import pytest
import warnings

from openpathsampling.deprecations import *

def make_foo_bar_baz():
    foo = Deprecation(problem="Foo is here.",
                      remedy="Get rid of the foo.",
                      remove_version=(2,0),
                      deprecated_in=(1,1,2))

    bar = Deprecation(problem="Bar is here.",
                      remedy="Get more foo!",
                      remove_version=(3,0),
                      deprecated_in=(0,5,4),
                      warn_once=False)

    baz = Deprecation(problem="Baz is too bazzerific.",
                      remedy="Find a qux.",
                      remove_version=(3,0),
                      deprecated_in=(1,3,2))

    return foo, bar, baz

FOO_DOCSTRING_MESSAGE = """

.. deprecated:: 1.1.2
    Foo is here. Get rid of the foo.
"""

def test_update_docstring():
    foo = make_foo_bar_baz()[0]
    def thing_with_docstring():
        """This is a docstring"""
        pass

    expected = "This is a docstring" + FOO_DOCSTRING_MESSAGE
    assert update_docstring(thing_with_docstring, foo) == expected

    def thing_without_docstring():
        pass

    assert update_docstring(thing_without_docstring, foo) == \
            FOO_DOCSTRING_MESSAGE

@pytest.mark.parametrize('version_tup,version_str',
                         [((2, 0), '2.0'), ((0, 9, 6), '0.9.6')],
                         ids=['2.0', '0.9.6'])
def test_version_tuple_to_string(version_tup, version_str):
    assert version_tuple_to_string(version_tup) == version_str

class TestDeprecation(object):
    def setup(self):
        self.foo, self.bar = make_foo_bar_baz()[0:2]

    def test_message(self):
        assert self.foo.message == "Foo is here. Get rid of the foo."

    def test_warn_once(self):
        with warnings.catch_warnings(record=True) as w:
            self.foo.warn()
            self.foo.warn()

        assert len(w) == 1
        assert issubclass(w[0].category, DeprecationWarning)
        assert self.foo.message in str(w[0].message)

    def test_warn_multiple(self):
        with warnings.catch_warnings(record=True) as w:
            self.bar.warn()
            self.bar.warn()

        assert len(w) == 2
        for warn in w:
            assert issubclass(warn.category, DeprecationWarning)
            assert self.bar.message in str(warn.message)


class TestDeprecationDecorators(object):
    def setup(self):
        self.foo, self.bar, self.baz = make_foo_bar_baz()

        @deprecate(self.foo)
        class RottedCode(object):
            __module__ = "fake_module"  # req'd for Py2

            @deprecate(self.bar)
            def bar(self):
                """Bar docstring"""
                pass

            # note the order of decorators here! important!
            @property
            @deprecate(self.baz)
            def baz(self):
                """Baz docstring"""
                pass

        self.RottedCode = RottedCode

    @pytest.mark.skipif(sys.version_info < (3,),
                        reason="Testing procedure requires Python 3")
    def test_docstrings(self):
        assert self.RottedCode.__doc__ is None
        assert self.RottedCode.bar.__doc__ == "Bar docstring"
        assert self.RottedCode.baz.__doc__ == "Baz docstring"

        self.RottedCode = has_deprecations(self.RottedCode)

        assert self.foo.message in self.RottedCode.__doc__
        assert self.bar.message in self.RottedCode.bar.__doc__
        # TODO: this does NOT work for properties (yet)
        # assert self.baz.message in self.RottedCode.baz.__doc__

    @pytest.mark.skipif(sys.version_info < (3,),
                        reason="Testing procedure requires Python 3")
    @pytest.mark.parametrize('use_has_deprecations', [True, False])
    def test_warnings(self, use_has_deprecations):
        wrapper = {True: has_deprecations,
                   False: lambda x: x}[use_has_deprecations]
        self.RottedCode = wrapper(self.RottedCode)

        with warnings.catch_warnings(record=True) as warns:
            rotted = self.RottedCode()
            rotted.bar()
            rotted.baz

        assert len(warns) == 3
        assert self.foo.message in str(warns[0].message)
        assert self.bar.message in str(warns[1].message)
        assert self.baz.message in str(warns[2].message)
        for w in warns:
            assert issubclass(w.category, DeprecationWarning)


def test_list_deprecations_default():
    # smoke test
    list_deprecations()

def test_list_deprecations_empty():
    assert list_deprecations(deprecations=[]) == []

def test_list_deprecations():
    deprecations = list(make_foo_bar_baz())
    foo, bar, baz = deprecations
    assert list_deprecations(deprecations=deprecations) == deprecations
    assert list_deprecations("2.5", deprecations) == [foo]
    assert list_deprecations("2.0", deprecations) == [foo]
    assert list_deprecations("4.0", deprecations) == [foo, bar, baz]
    assert list_deprecations("1.0", deprecations) == []
