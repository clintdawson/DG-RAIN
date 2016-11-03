import numpy as np
from math import sqrt

def read_channel_nodes(filename):
	"This function reads in the channel nodes and returns a list\
	containing the channels represented by their nodes."
	channelList = []
	dataDictionary = {}
	channelNodesFile = open(filename,"r")
	numChannels = int(channelNodesFile.readline())
	
	for i in range(0, numChannels):
		numNodes = int(channelNodesFile.readline())
		newChannel = []
		for j in range(0, numNodes):
			node, b, m1, m2 = [float(value) for value in channelNodesFile.readline().split()]
			node = int(node)
			node = node-1
			newChannel.append(node)
			#hold spot for nFriction (0.0 for now)
			dataArray = np.array([b, m1, m2,0.0])
			dataDictionary[node] = dataArray

		channelList.append(np.array(newChannel))
	
	channelNodesFile.close()
	return channelList, dataDictionary

def read_coordinates(filename, channelDictionary):
	r""" This function reads the coordinates of the nodes"""

	grid_file = open(filename, 'r')
	grid_file.readline()
	num_elements, num_nodes = [int(value) for value in grid_file.readline().split()]
	
	# Create necessary array storage
	coords =	{} 
	
	# Read in coordinates of each node
	for n in xrange(num_nodes):
	 	line = grid_file.readline().split()
		dataArray = np.array([float(line[1]), float(line[2]), float(line[3]), float(line[4])])
	 	coords[int(line[0])-1] = dataArray
	
	grid_file.close()

	# fill in the value of Manning's n
	for key in channelDictionary:
		channelAttrList = channelDictionary[key]
		channelAttrList[3] = coords[key][3]
		channelDictionary[key] = channelAttrList
	return coords


def find_junction_nodes(channelList):
	"This function finds the junction nodes"
	numChannels = len(channelList)
	junctionNodes = []

	# compare all channels with each other to find junction nodes
	for i in range(0,numChannels):
		for j in range(i+1,numChannels):
			arr1 = channelList[i]
			arr2 = channelList[j]
			commonNodes = np.intersect1d(arr1, arr2)
			for node in commonNodes:
				# add the node to the junction nodes list if it
				#doesn't already exist
				if node not in junctionNodes:
					junctionNodes.append(node)	
	return junctionNodes

def break_channels(junctionNodes, oldChannelList):
	"This function scans the channel nodes to ensure that a \
	junction doesn't occur in the middle of a channel. If it\
	does, then the function breaks the channels into smaller\
	channels"
	
	channelList = []
	junctionNodes = np.array(junctionNodes)
	# first make sure that a channel is connected to a junction
	# only at the end points
	for channel in oldChannelList:
		# find out how many junction nodes the channel contains
		indices = []
		intersection = np.intersect1d(junctionNodes, channel)
		if not intersection:
			channelList.append(channel)
		else:
			for node in intersection:
				indices.append(np.where(channel==node)[0][0])
			
			indices.sort()
			channelLength = len(channel)
			for i in range(0,len(indices)):
				if i == 0:
					firstIndex = 0;
				else:
					firstIndex = indices[i-1]
				
				subChannel1 = (channel[firstIndex:indices[i]+1]).tolist()
				subChannel2 = (channel[indices[i]:channelLength]).tolist()


				#subChannel1 = np.array(channel[firstIndex:indices[i]+1])
				#subChannel2 = np.array(channel[indices[i]:channelLength])
				if len(subChannel1) > 1:
					channelList.append(subChannel1)

				if len(subChannel2) > 1:
					channelList.append(subChannel2)

	return channelList




