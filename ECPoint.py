class ECPoint:
	def __init__(self,x, y = 0, z = 1):
		self.x = x
		self.y = y
		self.z = z

	def x(self):
		return self.x

	def z(self):
		return self.z

	def __str__(self):

		'''
		Converts a point P=(x,y) to P=(X,Z) 

		Params:
		P: point in affine coordinates

		Return:
		P in projectve coordinates
		'''
		
		return "({}:{})".format(self.x,self.z)

	def affine_coord(self):
		'''
		Converts a point P=(X,Z) to P=(x,y) 

		Params:
		P: point in projectve coordinates

		Return:
		P in affine coordinates
		'''

	# QUESTION: how do i find the right y coordinate? 

		# if P[1] != 0:
		# 	x = P[0] *  inverse(P[1], MOD)
		# 	y = 0xdeadbeef					# change this
		
		return "x(point) = " + (x * Util.mul_inverse(z))

	def inverse(self):
		pass