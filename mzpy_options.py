compression_methods = ["gamma", "golomb"]
compression_opts = ["rpos", "vpos", "poff"]

class CLIError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class CLIEmpty(Exception):
    pass
