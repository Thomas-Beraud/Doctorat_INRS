import numpy as np
import matplotlib.pyplot as plt

def gaussian(s, n, r, h):
	"""
	s : sill
	n : nugget
	r : range
	h : distance
	"""

	return (s-n)*(1-np.exp((-h**2)/(r**2*(1/3))))+n

def exponential(s, n, r, h):
	"""
	s : sill
	n : nugget
	r : range
	h : distance
	"""

	return (s-n)*(1-np.exp(-abs(h)/(r/3)))+n

def spherical(s, n, r, h):
	"""
	s : sill
	n : nugget
	r : range
	h : distance
	"""

	return (s-n)*(((3*abs(h))/(2*r)) - (abs(h)**3)/(2*(abs(h)**3)))+n



def clean_ellipse_index(list_of_index,value_of_var):
	"""
	This function take all the index around the ellipse (in a rectangular vicinity)
	and retain only the non-nill indexes to save computationnal time and memory
	"""

	non_zero_index = np.count_nonzero(value_of_var)

	temp_local_ellipse_index = [None] * non_zero_index
	temp_local_ellipse_value = np.zeros(non_zero_index,dtype=float)

	index = 0

	for item in range(len(list_of_index)):
		if value_of_var[item] != 0. :
			temp_local_ellipse_index[index] = list_of_index[item]
			temp_local_ellipse_value[index] = value_of_var[item]
			index += 1

	return(temp_local_ellipse_index, temp_local_ellipse_value)


def get_ellipse_bb(x, y, major, minor, angle_rad):
	"""

	major : vario_range along long axis
	minor : vario_range along short axis

	Compute tight ellipse bounding box.

	https://gist.github.com/smidm/b398312a13f60c24449a2c7533877dc0

	see https://stackoverflow.com/questions/87734/how-do-you-calculate-the-axis-aligned-bounding-box-of-an-ellipse#88020
	"""

	if angle_rad == 0 :

		return -major, -minor, major, minor

	else :

		t = np.arctan(-minor * np.tan(angle_rad) / (major))
		[max_x, min_x] = [x + major * np.cos(t) * np.cos(angle_rad) -
	  					minor * np.sin(t) * np.sin(angle_rad) for t in (t, t + np.pi)]

		t = np.arctan(minor * 1. / np.tan(angle_rad) / (major))
		[max_y, min_y] = [y + minor * np.sin(t) * np.cos(angle_rad) +
	  					major * np.cos(t) * np.sin(angle_rad) for t in (t, t + np.pi)]

		return np.floor(min_x), np.floor(min_y), np.ceil(max_x), np.ceil(max_y)



class Observation_Weight:
    def __init__(self,
                nx,
                ny,
                nz,
                x_step,
                y_step,
                z_step,
                max_range,
                med_range,
                min_range,
                azimuth,
                dip) :

        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.x_step = x_step
        self.y_step = y_step
        self.z_step = z_step
        self.max_range = max_range
        self.med_range = med_range
        self.min_range = min_range
        self.azimuth = azimuth
        self.dip = dip



        [self.x_min,self.y_min,self.x_max,self.y_max]=get_ellipse_bb(0,0, self.max_range,self.med_range,self.azimuth)
        #Compute 2D plan x-z
        [self.x_min,self.z_min,self.x_max,self.z_max]=get_ellipse_bb(0,0, self.max_range,self.min_range,self.dip)


        self.x_range = int(round(( self.x_max - self.x_min ) / self.x_step))+1 # Add one to count the 0 cell, symetric
        self.y_range = int(round(( self.y_max - self.y_min ) / self.y_step))+1 # Add one to count the 0 cell, symetric
        self.z_range = int(round(( self.z_max - self.z_min ) / self.z_step))+1 # Add one to count the 0 cell, symetric

        self.num_of_cells = int(self.x_range*self.y_range*self.z_range)

        self.local_ellipse_index = [None] * self.num_of_cells
        self.local_ellipse_value = np.zeros(self.num_of_cells,dtype=float)

        index = 0
        for z in range(self.z_range):
        	z_pos = z*self.z_step + self.z_min
        	for y in range(self.y_range):
        		y_pos = y*self.y_step + self.y_min
        		for x in range(self.x_range):
        			x_pos = x*self.x_step + self.x_min

        			self.local_ellipse_index[index] = [int(x+self.x_min/self.x_step),
                                                    int(y+self.y_min/self.y_step),
                                                    int(z+self.z_min/self.z_step)]


        			self.local_ellipse_value[index] = self.get_axis_range(x_pos,y_pos,z_pos,0,0,0,"exponential")
        			index += 1


        self.local_index, self.local_var = clean_ellipse_index(self.local_ellipse_index, self.local_ellipse_value)


    def get_range(self):

        return(self.z_range, self.y_range, self.x_range)

    """
    def position_to_index(self, x, y, z):


    	x, y, z : Coordinates of the observed cell in the model,
    				/ ! \ given in cell number and not real distance


    	return (x + y*self.nx_ + z*self.nx_*self.ny_)
        """


    def get_axis_range(self, x, y, z,x_c,y_c,z_c,variogram_type="exponential"):
    	"""
    	x_c, y_z and z_c : cell coordinates of the observed cell
    	x, y and c : coordinates of surroundings cells in the ellipsoid range search

    	variogram_type : String "gaussian", "exponential" or "spherical"
    	"""

    	s = 1.0
    	n = 0.0


    	coords_local = np.matrix([(x - x_c),(y - y_c),(z - z_c)])

    	# rotation_along_y to remove dipping

    	mat = np.matrix([[np.cos(self.dip), 0, np.sin(self.dip)],
    		 [0, 1, 0],
    		 [-np.sin(self.dip), 0, np.cos(self.dip)]])

    	coords_local = mat@coords_local.T

    	#rotation_along_z to remove azimuth

    	mat = np.matrix([[np.cos(self.azimuth), -np.sin(self.azimuth),0],
    	 [np.sin(self.azimuth), np.cos(self.azimuth), 0],
    	 [0, 0, 1]])

    	coords_local = mat@coords_local

    	# coordinates are now aligned with axis north, east, vertical axis without any angle_deg
    	# it is a necessity to compute a, b and c

    	a = np.abs(coords_local[0]/self.x_max)**2
    	b = np.abs(coords_local[1]/self.y_max)**2
    	c = np.abs(coords_local[2]/self.z_max)**2

    	if (a+b+c)<=1 :
    		if variogram_type == "exponential" :
    			return(s-(a*exponential(s,n,self.x_max,coords_local[0]) + \
                        b*exponential(s,n,self.y_max,coords_local[1]) + \
                        c*exponential(s,n,self.z_max,coords_local[2])))

    		if variogram_type == "gaussian" :
    			return(s-(a*gaussian(s,n,self.x_max,coords_local[0]) + \
                        b*gaussian(s,n,self.y_max,coords_local[1]) + \
                        c*gaussian(s,n,self.z_max,coords_local[2])))

    		if variogram_type == "spherical" :
    			return(s-(a*spherical(s,n,self.x_max,coords_local[0]) + \
                        b*spherical(s,n,self.y_max,coords_local[1]) + \
                        c*spherical(s,n,self.z_max,coords_local[2])))

    	if (a+b+c) == 0:
    		return(s-n)
    	else :
    		return(s-s)


    def global_ellispoid_research(self, local_index,local_var):
    	"""
    	matrix are size number of cells * 2 and maximum amount of non nil variance cell
    	"""

    	row = np.zeros((self.num_of_cells*len(local_var)),dtype=int)
    	col = np.zeros((self.num_of_cells*len(local_var)),dtype=int)
    	var = np.zeros((self.num_of_cells*len(local_var)),dtype=np.float32)

    	temp_index = 0
    	for index in range(self.num_of_cells):


    		z = index // (self.ny*self.nx_)
    		y = (index % (self.ny*self.nx) ) // self.nx
    		x = (index % (self.ny*self.nx) ) %  self.nx

    		for i in range(len(local_index)):
    			x_temp = x+local_index[i][0]
    			y_temp = y+local_index[i][1]
    			z_temp = z+local_index[i][2]

    			if x_temp >= 0 and x_temp < self.nx and y_temp >= 0 and y_temp < self.ny and z_temp >= 0 and z_temp < self.nz :
    				row[temp_index] = index
    				col[temp_index] = x_temp+y_temp*self.nx+z_temp*self.nx*self.ny
    				var[temp_index] = local_var[i]
    				temp_index += 1

    	end = np.count_nonzero(var)

    	return (row[:end],col[:end],var[:end])

    def ellispoid_around_index(self, x, y, z, local_index, local_var):

    	"""
    	index : Int giving the index of the cell in the model
    			/ ! \ Check the logic with your model, here number increases along x, y and finally z

    	local_index and local_var : non nill cell in the ellipsoid around the observed point,
    								obtained from the clean_ellipse_index function

    	nx, ny, nz : Dimensions of the model for a regular grid in grid number

    	return : Value of observation weight in all cells of the complete model

    	"""

    	var = np.zeros((self.num_of_cells,1),dtype=np.float32)

    	#z = index // (ny*nx)
    	#y = (index % (ny*nx) ) // nx
    	#x = (index % (ny*nx) ) %  nx

    	for i in range(len(local_index)):
    		x_temp = x+local_index[i][0]
    		y_temp = y+local_index[i][1]
    		z_temp = z+local_index[i][2]

    		if x_temp >= 0 and x_temp < self.nx and y_temp >= 0 and y_temp < self.ny and z_temp >= 0 and z_temp < self.nz :
    			temp_index = x_temp+y_temp*self.nx+z_temp*self.nx*self.ny
    			var[temp_index] = local_var[i]


    	return (var)
