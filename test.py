import numpy as np

def foo(a, b):
    # c is apparently pointing to the same memory that's holding the values in a
    # so when you assign c to a, then modify c, it actually modifies a !!
    c = a
    c[0] = 2

    # you can make a local variable called b, that is a copy of the content of the 
    # input variable b. then when you change some value in the local variable b, it has
    # no effect on the input variable b
    b = b.copy()
    b[0] = 2


print 'initialize '

a = np.zeros(10)
b = np.zeros(10)

print np.max(a)
print np.max(b)

foo(a,b)

print 'call function foo'

print np.max(a)
print np.max(b)


