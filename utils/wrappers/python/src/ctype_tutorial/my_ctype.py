// helpful tut: https://terryoy.github.io/2014/03/using-ctypes-to-wrap-c-library.html
// helpful tut: http://indico.ictp.it/event/7657/session/4/contribution/19/material/1/0.pdf
// ... the latter discusses passing arrays and structures

from ctypes import *
my_ctyp = CDLL("./lib_my_ctype.so")

a = c_int(3) # initialize a as a cytpe int with value 3
b = c_int(7)

my_ctyp.my_update(byref(a),byref(b)) # Pass them to my_update by ref

print a.value # Examine the new updated values
print b.value
