#!/usr/bin/env python3


class Master(object):
	def __init__(self):
		self.testval = 1

class Child(Master):
	def __init__(self):
		super(Child,self).__init__()
		self.hi = 'hi'
	def print_hi(self):
		print(self.hi)
		print(self.testval)


C = Child()
C.print_hi()

