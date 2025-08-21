
class HaPyException(Exception):
    """
    The spectral magnitude estimation class
    """

    def __init__(self, message="other"):
        self.message = "Mw estimation error: " + message
        super().__init__(self.message)


def print_separation_double_line():
    print('==================================================================')


def print_separation_single_line():
    print('------------------------------------------------------------------')
