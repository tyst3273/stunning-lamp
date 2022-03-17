
class a_class:

    def __init__(self,a=1):
        self.a = 1
        
    @property
    def a(self):
        return self._a

    @a.setter
    def a(self,value):
        print(value)
        self._a = value

    @a.deleter
    def a(self):
        print('deleting a')
        del self._a


if __name__ == "__main__":

    c = a_class()
    print(c.a)
    c.a = 'a'

    a = c.a 
    print(a)

    del c.a

