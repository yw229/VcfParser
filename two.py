#!/usr/bin/env python 

import one

print("top-level in two.py")

one.func()

if __name__ =="__main__":
	print("two.py is being called directly")
else:
	print("two.py is being imported into another module one")

