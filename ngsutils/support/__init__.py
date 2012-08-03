class Symbolize(object):
    'Converts strings to symbols - basically a cache of strings'
    def __init__(self):
        self.__cache = {}

    def __getitem__(self, k):
        if not k in self.__cache:
            self.__cache[k] = k

        return self.__cache[k]

symbols = Symbolize()
