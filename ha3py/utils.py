class HaPyException(Exception):
    """
    The spectral magnitude estimation exception class

    :param message: The exception message
    :type message: str

    """

    def __init__(self, message="other"):
        """
        :param message: The exception message
        :type message: str

        """
        self.message = "Mw estimation error: " + message
        super().__init__(self.message)


def print_separation_double_line():
    print('==================================================================')


def print_separation_single_line():
    print('------------------------------------------------------------------')
