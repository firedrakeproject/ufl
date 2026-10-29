"""Hashing of form arguments."""

from ufl import FunctionSpace, Mesh, TestFunction, interval

from utils import LagrangeElement


class DisplayFunctionSpace(FunctionSpace):
    """Equal spaces may carry different display-only details."""

    def __init__(self, domain, element, display):
        super().__init__(domain, element)
        self.display = display

    def __repr__(self):
        return f"DisplayFunctionSpace({self.display!r})"


def test_equal_arguments_have_equal_hashes_when_space_reprs_differ():
    domain = Mesh(LagrangeElement(interval, 1, (1,)))
    element = LagrangeElement(interval, 1)
    left = TestFunction(DisplayFunctionSpace(domain, element, "parent 0"))
    right = TestFunction(DisplayFunctionSpace(domain, element, "parent 1"))

    assert left == right
    assert hash(left) == hash(right)
    assert right in {left}
