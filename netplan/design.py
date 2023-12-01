#!/usr/bin/env python
"""
A heuristic approach for two-level network design - rural electrification
Contributors: Ayse Selin Kocaman ,ask2170@columbia.edu & Simone Fobi, sf2786@columbia.edu
Alex Zvoleff, aiz2101@columbia.edu
Ariel Nuñez, ariel@piensa.co
Justine Tunney, jart at github
The team at https://github.com/Turbo87/utm for utm functions (MIT License)
"""
import os
import csv
import sys
import time
import copy
import h3
import gc
import collections
import multiprocessing
import pathlib
import heapq
import math
import gzip
use_numpy = False

class OutOfRangeError(ValueError):
    pass

__all__ = ['to_latlon', 'from_latlon']

K0 = 0.9996

E = 0.00669438
E2 = E * E
E3 = E2 * E
E_P2 = E / (1 - E)

SQRT_E = math.sqrt(1 - E)
_E = (1 - SQRT_E) / (1 + SQRT_E)
_E2 = _E * _E
_E3 = _E2 * _E
_E4 = _E3 * _E
_E5 = _E4 * _E

M1 = (1 - E / 4 - 3 * E2 / 64 - 5 * E3 / 256)
M2 = (3 * E / 8 + 3 * E2 / 32 + 45 * E3 / 1024)
M3 = (15 * E2 / 256 + 45 * E3 / 1024)
M4 = (35 * E3 / 3072)

P2 = (3 / 2 * _E - 27 / 32 * _E3 + 269 / 512 * _E5)
P3 = (21 / 16 * _E2 - 55 / 32 * _E4)
P4 = (151 / 96 * _E3 - 417 / 128 * _E5)
P5 = (1097 / 512 * _E4)

R = 6378137

ZONE_LETTERS = "CDEFGHJKLMNPQRSTUVWXX"


def in_bounds(x, lower, upper, upper_strict=False):
    if upper_strict and use_numpy:
        return lower <= math.min(x) and math.max(x) < upper
    elif upper_strict and not use_numpy:
        return lower <= x < upper
    elif use_numpy:
        return lower <= math.min(x) and math.max(x) <= upper
    return lower <= x <= upper


def check_valid_zone(zone_number, zone_letter):
    if not 1 <= zone_number <= 60:
        raise OutOfRangeError('zone number out of range (must be between 1 and 60)')

    if zone_letter:
        zone_letter = zone_letter.upper()

        if not 'C' <= zone_letter <= 'X' or zone_letter in ['I', 'O']:
            raise OutOfRangeError('zone letter out of range (must be between C and X)')


def mixed_signs(x):
    return use_numpy and math.min(x) < 0 and math.max(x) >= 0


def negative(x):
    if use_numpy:
        return math.max(x) < 0
    return x < 0


def mod_angle(value):
    """Returns angle in radians to be between -pi and pi"""
    return (value + math.pi) % (2 * math.pi) - math.pi


def to_latlon(easting, northing, zone_number, zone_letter=None, northern=None, strict=True):
    """This function converts UTM coordinates to Latitude and Longitude

        Parameters
        ----------
        easting: int or NumPy array
            Easting value of UTM coordinates

        northing: int or NumPy array
            Northing value of UTM coordinates

        zone_number: int
            Zone number is represented with global map numbers of a UTM zone
            numbers map. For more information see utmzones [1]_

        zone_letter: str
            Zone letter can be represented as string values.  UTM zone
            designators can be seen in [1]_

        northern: bool
            You can set True or False to set this parameter. Default is None

        strict: bool
            Raise an OutOfRangeError if outside of bounds

        Returns
        -------
        latitude: float or NumPy array
            Latitude between 80 deg S and 84 deg N, e.g. (-80.0 to 84.0)

        longitude: float or NumPy array
            Longitude between 180 deg W and 180 deg E, e.g. (-180.0 to 180.0).


       .. _[1]: http://www.jaworski.ca/utmzones.htm

    """
    if not zone_letter and northern is None:
        raise ValueError('either zone_letter or northern needs to be set')

    elif zone_letter and northern is not None:
        raise ValueError('set either zone_letter or northern, but not both')

    if strict:
        if not in_bounds(easting, 100000, 1000000, upper_strict=True):
            raise OutOfRangeError('easting out of range (must be between 100,000 m and 999,999 m)')
        if not in_bounds(northing, 0, 10000000):
            raise OutOfRangeError('northing out of range (must be between 0 m and 10,000,000 m)')
    
    check_valid_zone(zone_number, zone_letter)
    
    if zone_letter:
        zone_letter = zone_letter.upper()
        northern = (zone_letter >= 'N')

    x = easting - 500000
    y = northing

    if not northern:
        y -= 10000000

    m = y / K0
    mu = m / (R * M1)

    p_rad = (mu +
             P2 * math.sin(2 * mu) +
             P3 * math.sin(4 * mu) +
             P4 * math.sin(6 * mu) +
             P5 * math.sin(8 * mu))

    p_sin = math.sin(p_rad)
    p_sin2 = p_sin * p_sin

    p_cos = math.cos(p_rad)

    p_tan = p_sin / p_cos
    p_tan2 = p_tan * p_tan
    p_tan4 = p_tan2 * p_tan2

    ep_sin = 1 - E * p_sin2
    ep_sin_sqrt = math.sqrt(1 - E * p_sin2)

    n = R / ep_sin_sqrt
    r = (1 - E) / ep_sin

    c = E_P2 * p_cos**2
    c2 = c * c

    d = x / (n * K0)
    d2 = d * d
    d3 = d2 * d
    d4 = d3 * d
    d5 = d4 * d
    d6 = d5 * d

    latitude = (p_rad - (p_tan / r) *
                (d2 / 2 -
                 d4 / 24 * (5 + 3 * p_tan2 + 10 * c - 4 * c2 - 9 * E_P2)) +
                 d6 / 720 * (61 + 90 * p_tan2 + 298 * c + 45 * p_tan4 - 252 * E_P2 - 3 * c2))

    longitude = (d -
                 d3 / 6 * (1 + 2 * p_tan2 + c) +
                 d5 / 120 * (5 - 2 * c + 28 * p_tan2 - 3 * c2 + 8 * E_P2 + 24 * p_tan4)) / p_cos

    longitude = mod_angle(longitude + math.radians(zone_number_to_central_longitude(zone_number)))

    return (math.degrees(latitude),
            math.degrees(longitude))


def from_latlon(latitude, longitude, force_zone_number=None, force_zone_letter=None):
    """This function converts Latitude and Longitude to UTM coordinate

        Parameters
        ----------
        latitude: float or NumPy array
            Latitude between 80 deg S and 84 deg N, e.g. (-80.0 to 84.0)

        longitude: float or NumPy array
            Longitude between 180 deg W and 180 deg E, e.g. (-180.0 to 180.0).

        force_zone_number: int
            Zone number is represented by global map numbers of an UTM zone
            numbers map. You may force conversion to be included within one
            UTM zone number.  For more information see utmzones [1]_

        force_zone_letter: str
            You may force conversion to be included within one UTM zone
            letter.  For more information see utmzones [1]_

        Returns
        -------
        easting: float or NumPy array
            Easting value of UTM coordinates

        northing: float or NumPy array
            Northing value of UTM coordinates

        zone_number: int
            Zone number is represented by global map numbers of a UTM zone
            numbers map. More information see utmzones [1]_

        zone_letter: str
            Zone letter is represented by a string value. UTM zone designators
            can be accessed in [1]_


       .. _[1]: http://www.jaworski.ca/utmzones.htm
    """
    if not in_bounds(latitude, -80, 84):
        raise OutOfRangeError('latitude out of range (must be between 80 deg S and 84 deg N)')
    if not in_bounds(longitude, -180, 180):
        raise OutOfRangeError('longitude out of range (must be between 180 deg W and 180 deg E)')
    if force_zone_number is not None:
        check_valid_zone(force_zone_number, force_zone_letter)

    lat_rad = math.radians(latitude)
    lat_sin = math.sin(lat_rad)
    lat_cos = math.cos(lat_rad)

    lat_tan = lat_sin / lat_cos
    lat_tan2 = lat_tan * lat_tan
    lat_tan4 = lat_tan2 * lat_tan2

    if force_zone_number is None:
        zone_number = latlon_to_zone_number(latitude, longitude)
    else:
        zone_number = force_zone_number

    if force_zone_letter is None:
        zone_letter = latitude_to_zone_letter(latitude)
    else:
        zone_letter = force_zone_letter

    lon_rad = math.radians(longitude)
    central_lon = zone_number_to_central_longitude(zone_number)
    central_lon_rad = math.radians(central_lon)

    n = R / math.sqrt(1 - E * lat_sin**2)
    c = E_P2 * lat_cos**2

    a = lat_cos * mod_angle(lon_rad - central_lon_rad)
    a2 = a * a
    a3 = a2 * a
    a4 = a3 * a
    a5 = a4 * a
    a6 = a5 * a

    m = R * (M1 * lat_rad -
             M2 * math.sin(2 * lat_rad) +
             M3 * math.sin(4 * lat_rad) -
             M4 * math.sin(6 * lat_rad))

    easting = K0 * n * (a +
                        a3 / 6 * (1 - lat_tan2 + c) +
                        a5 / 120 * (5 - 18 * lat_tan2 + lat_tan4 + 72 * c - 58 * E_P2)) + 500000

    northing = K0 * (m + n * lat_tan * (a2 / 2 +
                                        a4 / 24 * (5 - lat_tan2 + 9 * c + 4 * c**2) +
                                        a6 / 720 * (61 - 58 * lat_tan2 + lat_tan4 + 600 * c - 330 * E_P2)))

    if mixed_signs(latitude):
        raise ValueError("latitudes must all have the same sign")
    elif negative(latitude):
        northing += 10000000

    return easting, northing, zone_number, zone_letter


def latitude_to_zone_letter(latitude):
    # If the input is a numpy array, just use the first element
    # User responsibility to make sure that all points are in one zone
    if use_numpy and isinstance(latitude, math.ndarray):
        latitude = latitude.flat[0]

    if -80 <= latitude <= 84:
        return ZONE_LETTERS[int(latitude + 80) >> 3]
    else:
        return None


def latlon_to_zone_number(latitude, longitude):
    # If the input is a numpy array, just use the first element
    # User responsibility to make sure that all points are in one zone
    if use_numpy:
        if isinstance(latitude, math.ndarray):
            latitude = latitude.flat[0]
        if isinstance(longitude, math.ndarray):
            longitude = longitude.flat[0]

    if 56 <= latitude < 64 and 3 <= longitude < 12:
        return 32

    if 72 <= latitude <= 84 and longitude >= 0:
        if longitude < 9:
            return 31
        elif longitude < 21:
            return 33
        elif longitude < 33:
            return 35
        elif longitude < 42:
            return 37

    return int((longitude + 180) / 6) + 1


def zone_number_to_central_longitude(zone_number):
    return (zone_number - 1) * 6 - 180 + 3


class Network():
    'Implements an undirected network made up of segments and nodes.'
    def __init__(self):
        self._network = {}
        self._netIDByNode = {}
        self._nodesByNetID = {}

    def __repr__(self):
        #TODO: Finish this
        return "__repr__ UNDEFINED"
        
    def __str__(self):
        #TODO: Finish this
        return "__str__ UNDEFINED"

    def addSeg(self, newSeg):
        'Adds a new segment to a network, taking care of netID assignment.'
        node1 = newSeg.getNode1()
        node2 = newSeg.getNode2()
        # Neither node already in network:
        if (not self.inNet(node1)) and (not self.inNet(node2)):
            try:
                newID = max(self.listNetIDs()) + 1 # fails for first seg
            except ValueError:
                # only occurs when first two nodes are added to the network
                newID = 0;
            self._netIDByNode[node1] = newID
            self._netIDByNode[node2] = newID
            self._nodesByNetID[newID] = [node1, node2]
            self._network[newID] = [newSeg]
        # Both or only one node already in network:
        else:
            # Both nodes already in network:
            if self.inNet(node1) and self.inNet(node2):
                baseNetID = min(self.getNetID(node1), self.getNetID(node2))
                mergingNetID = max(self.getNetID(node1), self.getNetID(node2))
                mergingNodes = self._nodesByNetID[mergingNetID]
                for node in mergingNodes:
                    self._netIDByNode[node] = baseNetID
                self._network[baseNetID].extend(
                        self._network.pop(mergingNetID))
                self._nodesByNetID[baseNetID].extend(
                        self._nodesByNetID.pop(mergingNetID))
            # Only one node already in network:
            else:
                if self.inNet(node1):
                    baseNode = node1
                    newNode = node2
                else:
                    baseNode = node2
                    newNode = node1
                baseNetID = self.getNetID(baseNode)
                self._nodesByNetID[baseNetID].append(newNode)
                self._netIDByNode[newNode] = baseNetID
            self._network[baseNetID].append(newSeg)
        return 0

    def inNet(self, value):
        if value in self._netIDByNode.keys():
            return True
        else:
            return False

    # def setInNetValue(self, node, value):
        # self._inNet(node)=value

    def getNetID(self, node):
        if node in self._netIDByNode.keys():
            return self._netIDByNode[node]
        else:
            raise ValueError("Node not in network")
        
    def getEdges(self):
        netIDs = self.listNetIDs()
        edgeList = []
        for netID in netIDs:
            for edge in self._network[netID]:
                edgeList.append(edge)
        return edgeList
    
    def getNodes(self):
        return self._netIDByNode.keys()

    def numNodes(self):
        return len(self._netIDByNode.keys())

    def listNetIDs(self):
        return self._nodesByNetID.keys()

    def getTotalEdgeWeight(self):
        totalWeight = 0
        for edge in self.getEdges():
            totalWeight += edge.getWeight()
        return totalWeight


class Seg:
    'A class representing undirected segs.'
    def __init__(self, ID=None, node1=None, node2=None, weight=None):
        self._ID = ID
        self._node1 = node1
        self._node2 = node2
        self._weight = weight

    def __eq__(self, other):
        return self._ID == other._ID

    def __ne__(self, other):
        return self._ID != other._ID

    def __lt__(self, other):
        return self._weight < other._weight

    def __le__(self, other):
        return self._weight <= other._weight

    def __gt__(self, other):
        return self._weight > other._weight

    def __ge__(self, other):
        return self._weight >= other._weight

    def __hash__(self):
        return hash(self._ID)


    def getWeight(self):
        return self._weight
    
    def getID(self):
        return self._ID

    def getNode1(self):
        return self._node1

    def getNode2(self):
        return self._node2

    def getNodes(self):
        return self._node1, self._node2
    
    
class Node:
    'Defines a node class, with ID, x, y, weight and demand attributes'
    def __init__(self, value1=None, value2=None, value3=None, value4=None):
        self._id = value1
        self._x = value2
        self._y = value3
        self._weight = value4

    def __hash__(self):
        return hash(self._id)

    def __eq__(self, other):
        return self._id == other._id

    def __ne__(self, other):
        return self._id != other._id

    def __lt__(self, other):
        return self._weight < other._weight

    def __le__(self, other):
        return self._weight <= other._weight

    def __gt__(self, other):
        return self._weight > other._weight

    def __ge__(self, other):
        return self._weight >= other._weight

    def getID(self):
        return self._id

    def setID(self,newID):
        self._id=newID

    def getX(self):
        return self._x

    def setXY(self,X,Y):
        self._x=X
        self._y=Y

    def getY(self):
        return self._y

    def getWeight(self):
        return self._weight
    
    def setWeight(self,newWeight): # newWeight is NewMVMax or new distance in CMST
        self._weight=newWeight


class T:
    def __init__(self,value1=None,value2=None,value3=None,value4=None):
        self._removedSegmentNode= value1
        self._addedSegmentNode1=value2
        self._addedSegmentNode2=value3
        self._value=value4


    def __str__(self):
        value=self.getValue()
        return "T ( %(removedSegment)s), %(addedSegmentNode1)s, %(addedSegmentNode2)s, %(value)s)" %vars()

    def getRemovedSegmentNode(self):
        return self._removedSegmentNode
    def getAddedSegmentNode1(self):
        return self._addedSegmentNode1
    def getAddedSegmentNode2(self):
        return self._addedSegmentNode2
    def getValue(self):
        return self._value

def buildAssocDict(segments):
    'Builds dictionary with nodeID key where values are segs from/to that node'
    segList = {}
    for seg in segments:
        node1, node2 = seg.getNodes()
        for nodeID in [node1.getID(), node2.getID()]:
            #if segList.has_key(nodeID):
            if nodeID in segList.keys():
                segList[nodeID].append(seg)
            else:
                segList[nodeID] = [seg]
    return segList    


def depthFirstSet(node1,node2,root,segList,distDict):
    """
    Given a starting vertex, root, do a depth-first search and set the weights.
    """
    to_visit = []  # a list can be used as a stack in Python
    visited=[]
    node1.setWeight(node2.getWeight()+distDict[(node1.getID(),node2.getID())])
    to_visit.append(node1) # Start with root
    while len(to_visit)!= 0:
        v = to_visit.pop()
        #print v
        if v not in visited:
            visited.append(v)
            vNeighbors=[]
            #print segList[v.getID()]
            for seg in segList[v.getID()]:
                firstNode,secondNode=seg.getNodes()
                if firstNode==root or secondNode==root :
                    continue
                if firstNode==v and secondNode not in visited:
                    vNeighbors.append(secondNode)
                    secondNode.setWeight(v.getWeight()+distDict[(v.getID(),secondNode.getID())])
                if secondNode==v and firstNode not in visited:
                    vNeighbors.append(firstNode)
                    firstNode.setWeight(v.getWeight()+distDict[(v.getID(),firstNode.getID())])
            to_visit.extend(vNeighbors)

def maxDistFromNode(node1,root,treeSegments,distDict,segList):
    """
    Given a starting vertex, do a depth-first search and find the furthest node to that vertex
    """
    to_visit = []  # a list can be used as a stack in Python
    visited = set()

    tempWeightByNode={}
    to_visit.append(node1)
    tempWeightByNode[node1.getID()]=0
    while to_visit:
        v = to_visit.pop()
        if v not in visited:
            visited.add(v)
            vNeighbors=[]
            for seg in segList[v.getID()]:
                firstNode,secondNode=seg.getNodes()
                if firstNode==root or secondNode==root :
                    continue
                if firstNode==v and secondNode not in visited:
                    vNeighbors.append(secondNode)
                    tempWeightByNode[secondNode.getID()]=tempWeightByNode[v.getID()]+distDict[(v.getID(),secondNode.getID())]
                if secondNode==v and firstNode not in visited:
                    vNeighbors.append(firstNode)
                    tempWeightByNode[firstNode.getID()]=tempWeightByNode[v.getID()]+distDict[(v.getID(),firstNode.getID())]
            to_visit.extend(vNeighbors)
    maxDist=max(tempWeightByNode.values())
    return maxDist


            
def CMST(households,capacity,root):
    #Connects hhs directly to the root first
    treeSegments={}
    distDict={}
    households_Copy=copy.deepcopy(households)
    SegID=10000000
    root_Copy=copy.deepcopy(root)
    newRootID=root_Copy.getID()*(-1)-100  #### not to be confused with the same nodeID
    root_Copy.setID(newRootID)
    maxTvalue=0
    branchNodeByNode={}# which branch is the node on?
    nodesByBranchNode=collections.defaultdict(list )# what are the nodes on a spesific branch?
    for node in households_Copy:
        length=((node.getX()-root_Copy.getX())**2+(node.getY()-root_Copy.getY())**2)**(.5)
        treeSegments[(node.getID(),newRootID)]=Seg(SegID, node, root_Copy, length)
        SegID+=1
        node.setWeight(length)
        branchNodeByNode[node]=node
        nodesByBranchNode[node].append(node)
        distDict[(newRootID,node.getID())]=length
        distDict[(node.getID(),newRootID)]=length
       
    for node1 in households_Copy:
        for node2 in households_Copy:
            if node1==node2:
                continue
            else:
                distance=((node1.getX()-node2.getX())**2+(node1.getY()-node2.getY())**2)**(.5)
                distDict[(node1.getID(), node2.getID())]=distance
                Tvalue=treeSegments[(node1.getID(),newRootID)].getWeight()-distance

                if ((node2.getWeight()+distance)<=capacity):
                    newT=T(node1,node1,node2,Tvalue)

                    
                    if (newT.getValue()>=maxTvalue):
                        maxTObject=newT
                        maxTvalue=newT.getValue()
    totalLVCost=0
    for segment in treeSegments.values():
        totalLVCost=totalLVCost+segment.getWeight() 
    while(maxTvalue>0):

        maxTvalue=0
        SegID+=1
        node1=maxTObject.getAddedSegmentNode1() # node1 and node2 of new segment
        #print "node1", node1
        node2=maxTObject.getAddedSegmentNode2()
        #print "node2", node2
        
        #delete root segment of branch
        del treeSegments[(maxTObject.getRemovedSegmentNode().getID(),newRootID)]
       
        
        #if node1 is the first node in the branch
        if node1==branchNodeByNode[node1]:
            tempWeight=node1.getWeight() # I need this becouse node1 is updated first and it effects the others
            for node in nodesByBranchNode[branchNodeByNode[node1]]:
                node.setWeight(node.getWeight()-tempWeight+ node2.getWeight()+distDict[(node1.getID(),node2.getID())])# digerlerinin de set olmasi lazim

        #if not things get complicated and we need dfs
        else:
            segList=buildAssocDict(treeSegments.values()) # daha efficient yapmak icin sadece update edebilirim
            depthFirstSet(node1,node2,root_Copy,segList,distDict) # root  (node1) hala icinde unutma


        #tree updated after weights are set
        treeSegments[(node1.getID(), node2.getID())]=Seg(SegID,node1, node2,distDict[(node1.getID(),node2.getID())])
        segList = buildAssocDict(treeSegments.values())


        # Update dictionaries
        nodesByBranchNode[branchNodeByNode[node2]].extend(nodesByBranchNode.pop(branchNodeByNode[node1]))
        for node in nodesByBranchNode[branchNodeByNode[node2]]:
           branchNodeByNode[node]=branchNodeByNode[node2]

        # Rebuilt TStore & select maxT object
        for node1 in households_Copy:  #
            #print "node1", node1
            maxDistFromNode1=maxDistFromNode(node1,root_Copy,treeSegments,distDict,segList)
            for node2 in households_Copy:
                if node1==node2 or branchNodeByNode[node1]==branchNodeByNode[node2]:
                    continue
                else:
                    #print "maaxx",maxDistFromNode1
                    if (node2.getWeight()+distDict[node1.getID(),node2.getID()]+maxDistFromNode1<=capacity): #1 2ye baslansa ne olur?
                        #print "TTTTT", Tvalue
                        Tvalue=treeSegments[(branchNodeByNode[node1].getID(),newRootID)].getWeight()-distDict[(node1.getID(),node2.getID())]
                        if Tvalue>=maxTvalue:
                            maxTvalue=Tvalue
                            maxTObject=T(branchNodeByNode[node1],node1,node2,Tvalue)
        #print maxTvalue
        #print maxTObject
    totalLVCost=0    
    for segment in treeSegments.values():
        totalLVCost=totalLVCost+segment.getWeight()
    del root_Copy
    gc.collect()
    #print treeSegments.keys()
    #print households_Copy
    return treeSegments, totalLVCost


def generateSegments(centers, searchRadius):
    segments = []
    nodeCopy = centers.copy()

    segID = 0
    for startNode in centers.values():
        del nodeCopy[startNode.getID()]
        for endNode in nodeCopy.values():
            dist = ((startNode.getX() - endNode.getX()) ** 2 +
                    (startNode.getY() - endNode.getY()) ** 2) ** (.5)
            if dist < searchRadius:
                segments.append(Seg(segID, startNode, endNode, dist))
                segID += 1
    return segments


def maxTempInClusterDist(segment, ClusterByNode, nodesByClusterID):
    maxDist = 0

    tempCenter1, tempCenter2 = segment.getNodes()

    tempCenterX = (tempCenter1.getWeight() * tempCenter1.getX()
                   + tempCenter2.getWeight() * tempCenter2.getX()) / (tempCenter2.getWeight() + tempCenter1.getWeight())
    tempCenterY = (tempCenter1.getWeight() * tempCenter1.getY()
                   + tempCenter2.getWeight() * tempCenter2.getY()) / (tempCenter2.getWeight() + tempCenter1.getWeight())
    for node in nodesByClusterID[ClusterByNode[segment.getNode1()]]:
        dist = ((tempCenterX - node.getX()) ** 2 + (tempCenterY - node.getY()) ** 2) ** (.5)
        if dist >= maxDist:
            maxDist = dist

    for node in nodesByClusterID[ClusterByNode[segment.getNode2()]]:
        dist = ((tempCenterX - node.getX()) ** 2 + (tempCenterY - node.getY()) ** 2) ** (.5)
        if dist >= maxDist:
            maxDist = dist

    return maxDist, tempCenterX, tempCenterY


def loggers(log_filename, initial=False, data=None):
    if initial:
        with open(log_filename, 'w') as src:
            src.write("")
    else:
        with open(log_filename, 'a') as src:
            src.write(str(data) + "\n")



def primsAlg(segments, numNodes, firstNodeID, nodeDict):
    'Prim\'s Algorithm for finding a minimum spanning tree'

    tree = Network()
    segHeap = []

    # Find the shortest segment emanating from the node with the firstNodeID
    try:
        segs = nodeDict[firstNodeID]
    except KeyError:
        return tree

    leastWeight = None
    for seg in segs:
        if (leastWeight == None):
            leastWeight = seg.getWeight()
            firstSeg = seg
        elif (seg.getWeight() < leastWeight):
            leastWeight = seg.getWeight()
            firstSeg = seg
    tree.addSeg(firstSeg)

    # Starter to algorithm
    # Add the segs emanating from the first two endpoints to the heap
    for endNode in [firstSeg.getNode1(), firstSeg.getNode2()]:
        addToHeap(segHeap, nodeDict[endNode.getID()])

    # Pick best from heap and repeat
    while tree.numNodes() < numNodes:
        try:
            # Get best segment from heap
            seg = heapq.heappop(segHeap)
        except:
            # Tree is finished (not all nodes contained).
            break
        node1, node2 = seg.getNodes()
        node1InNet = tree.inNet(node1)
        node2InNet = tree.inNet(node2)
        # Add the seg if it's terminal node isn't already in the cluster.
        if (not node1InNet) or (not node2InNet):
            if not node1InNet:
                endNode = node1
            else:
                endNode = node2
            tree.addSeg(seg)
            # Add all emanating segs to the heap:
            # nodeDict returns all segments coming out from the endNode
            # endNode is the node that is outside of the tree
            addToHeap(segHeap, nodeDict[endNode.getID()])
            # And we are sure that everything in the heap is adjacent to the tree because
            # we only add the adjacent segments in the first place using nodeDict
    return tree


def addToHeap(heap, newSegs):
    'Adds new segments to the segHeap.'
    for seg in newSegs:
        heapq.heappush(heap, seg)
    return heap


def buildAssocDict(segments):
    'Builds dictionary with nodeID key where values are segs from/to that node'
    segList = {}
    for seg in segments:
        node1, node2 = seg.getNodes()
        for nodeID in [node1.getID(), node2.getID()]:
            if nodeID in segList.keys():
            # if segList.has_key(nodeID):
                segList[nodeID].append(seg)
            else:
                segList[nodeID] = [seg]
    return segList


def run(centers, nodesByClusterID, clusterByNode, LVCostDict, sr, MV, LV, TCost, distFromT, maxLVLenghtInCluster,
        outputDir, logfilename):
    sumLVCostAtEachStep = {}
    minCenters = copy.deepcopy(centers)
    segments = generateSegments(minCenters, sr)

    # To write total cost to a text file
    statFile = outputDir + os.sep + "TotalCost_FirstStage.txt"
    outFile = open(statFile, "w")

    minTotalCost = len(centers) * TCost
    outFile.write("%(minTotalCost)f\n" % vars())
    minLVCostDict = copy.deepcopy(LVCostDict)

    minNodesByClusterID = copy.deepcopy(nodesByClusterID)
    minClusterByNode = copy.deepcopy(clusterByNode)
    minSeg = min(segments, key=lambda obj: obj.getWeight())

    if minSeg.getWeight() <= distFromT * 2:
        maxDist = 0
    else:
        maxDist = distFromT + 10
        #print("NO CLUSTER POSSIBLE")

    tempCenter1, tempCenter2 = minSeg.getNodes()
    tempCenterX = (tempCenter1.getWeight() * tempCenter1.getX()
                   + tempCenter2.getWeight() * tempCenter2.getX()) / (tempCenter2.getWeight() + tempCenter1.getWeight())
    tempCenterY = (tempCenter1.getWeight() * tempCenter1.getY()
                   + tempCenter2.getWeight() * tempCenter2.getY()) / (tempCenter2.getWeight() + tempCenter1.getWeight())
    i = len(centers)
    initial = True
    loggers(logfilename, initial)
    initial = False
    while (maxDist <= distFromT):
        i -= 1
        cur_token = 'stage1 ' + str(i)
        loggers(logfilename, initial, cur_token)

        center1, center2 = minSeg.getNodes()
        weight = center2.getWeight() + center1.getWeight()
        baseClusterID = min(clusterByNode[center1], clusterByNode[center2])
        mergingClusterID = max(clusterByNode[center1], clusterByNode[center2])
        nodesByClusterID[baseClusterID].extend(nodesByClusterID.pop(mergingClusterID))

        centers[baseClusterID].setXY(tempCenterX, tempCenterY)
        centers[baseClusterID].setWeight(weight)

        del centers[mergingClusterID]

        for node in nodesByClusterID[baseClusterID]:
            clusterByNode[node] = baseClusterID

        segments = generateSegments(centers, sr)
        TotalTransformerCost = len(centers) * TCost

        del LVCostDict[mergingClusterID]
        gc.collect()
        _, LVCostDict[baseClusterID] = CMST(nodesByClusterID[baseClusterID],
                                                                    maxLVLenghtInCluster,
                                                                    centers[baseClusterID])
        # sums the cost
        sumLVCostAtEachStep[len(centers)] = sum(LVCostDict.values()) * LV
        newTotalCost = TotalTransformerCost + (sum(LVCostDict.values())) * LV

        outFile.write("%i %f\n" % (i, sumLVCostAtEachStep[len(centers)]))
        if (newTotalCost <= minTotalCost):
            minNodesByClusterID = copy.deepcopy(nodesByClusterID)
            # minTree=copy.deepcopy(newTree)
            minCenters = copy.deepcopy(centers)
            minLVCostDict = LVCostDict.copy()
            minTotalCost = newTotalCost
            minClusterByNode = copy.deepcopy(clusterByNode)
        # Calculate maxDist below for next graph and continue if it is less than 500

        try:  # to check if there is a segment on the graph or there is only one cluster  # bir tane break eden varsa bile devamini check ediyor!!!!!
            # seems this looks for the shortest segment with the lv less that distFromT
            minSeg = min(segments) # finds closest 2 points
            maxDist, tempCenterX, tempCenterY = maxTempInClusterDist(minSeg, clusterByNode, nodesByClusterID) # find the largest distance to centroid of minSeg
            if maxDist > distFromT: # if largest distance > 500 meters
                segments.sort() # sort by smallest to largest distances

                for seg in segments:
                    if seg.getWeight() > distFromT * 2: # distance greater than 1000 skip
                        break
                    else: # if distance is okay check if there are 2 closer points
                        maxDist, tempCenterX, tempCenterY = maxTempInClusterDist(seg, clusterByNode, nodesByClusterID)
                        if maxDist <= distFromT:
                            minSeg = seg  ## identifies a new minSeg to go to if there is still room to add to the LV
                            break # finds the next minimum segment to go to
        except:
            break

    outFile.close()

    if len(minCenters) == len(centers) or len(minCenters) == 1:
        segments_ST = generateSegments(centers, sr)
        nodeDict = buildAssocDict(segments_ST)
        minTree = primsAlg(segments_ST, len(centers), 0, nodeDict)
        minTotalCost_ST = minTree.getTotalEdgeWeight() * MV + minTotalCost
        return minTotalCost_ST, minTree, centers, nodesByClusterID, sum(LVCostDict.values()) * LV

    centers_ST = copy.deepcopy(minCenters)
    LVCostDict_ST = LVCostDict


    segments_ST = generateSegments(centers_ST, sr)

    # To write total cost to a text file
    statFile = outputDir + os.sep + "TotalCost_SecondStage.txt"
    outFile = open(statFile, "w")

    nodeDict = buildAssocDict(segments_ST)

    minTree = primsAlg(segments_ST, len(centers_ST), 0, nodeDict)  # 0 is the starting node of Prims algorithm
    i = len(centers_ST)
    minTotalCost_ST = minTree.getTotalEdgeWeight() * MV + len(centers_ST) * TCost + (sum(minLVCostDict.values())) * LV
    outFile.write(
        "%i %f %f %f\n" % (i, (sum(minLVCostDict.values())) * LV, minTree.getTotalEdgeWeight() * MV, minTotalCost_ST))

    minLVCostSum_ST = 9999999999999999  # a big number
    nodesByClusterID_ST = copy.deepcopy(minNodesByClusterID)
    clusterByNode_ST = copy.deepcopy(minClusterByNode)
    try:
        # given the tx location
        minSeg_ST = min(segments_ST, key=lambda obj: obj.getWeight()) #
        maxDist, tempCenterX, tempCenterY = maxTempInClusterDist(minSeg_ST, clusterByNode_ST, nodesByClusterID_ST)
        if maxDist > distFromT:
            segments_ST.sort(key=lambda obj: obj.getWeight())

            for seg in segments_ST:
                if seg.getWeight() > distFromT * 2:
                    break  # break from for loop
                else:
                    maxDist, tempCenterX, tempCenterY = maxTempInClusterDist(seg, clusterByNode_ST, nodesByClusterID_ST)
                    if maxDist <= distFromT:
                        minSeg_ST = seg
                        break  # break from for loop
    except:
        return minTotalCost_ST, minTree, centers_ST, nodesByClusterID_ST, sum(LVCostDict_ST.values()) * LV

    minNodesByClusterID_ST = copy.deepcopy(minNodesByClusterID)
    minCenters_ST = copy.deepcopy(minCenters)

    if minSeg_ST.getWeight() <= distFromT * 2:
        maxDist = 0
    else:
        maxDist = distFromT + 10
        print("NO CLUSTER POSSIBLE")

    tempCenter1, tempCenter2 = minSeg_ST.getNodes()

    tempCenterX = (tempCenter1.getWeight() * tempCenter1.getX()
                   + tempCenter2.getWeight() * tempCenter2.getX()) / (tempCenter2.getWeight() + tempCenter1.getWeight())
    tempCenterY = (tempCenter1.getWeight() * tempCenter1.getY()
                   + tempCenter2.getWeight() * tempCenter2.getY()) / (tempCenter2.getWeight() + tempCenter1.getWeight())

    initial = False
    i = len(minCenters)
    while (maxDist <= distFromT):

        i -= 1
        if i % 20 == 0:
            cur_token = 'stage2 ' + str(i)
            loggers(logfilename, initial, cur_token)
        center1, center2 = minSeg_ST.getNodes()

        weight = center2.getWeight() + center1.getWeight()
        baseClusterID = min(clusterByNode_ST[center1], clusterByNode_ST[center2])

        mergingClusterID = max(clusterByNode_ST[center1], clusterByNode_ST[center2])

        nodesByClusterID_ST[baseClusterID].extend(nodesByClusterID_ST.pop(mergingClusterID))

        centers_ST[baseClusterID].setXY(tempCenterX, tempCenterY)
        centers_ST[baseClusterID].setWeight(weight)

        del centers_ST[mergingClusterID]

        for node in nodesByClusterID_ST[baseClusterID]:
            clusterByNode_ST[node] = baseClusterID

        segments_ST = generateSegments(centers_ST, sr)
        nodeDict = buildAssocDict(segments_ST)
        newTree = primsAlg(segments_ST, len(centers_ST), 0, nodeDict)
        TotalMVCost_ST = newTree.getTotalEdgeWeight() * MV
        TotalTransformerCost_ST = len(centers_ST) * TCost
        gc.collect()

        newTotalCost_ST = TotalMVCost_ST + TotalTransformerCost_ST + sumLVCostAtEachStep[len(centers_ST)]

        if (newTotalCost_ST <= minTotalCost_ST):
            minNodesByClusterID_ST = copy.deepcopy(nodesByClusterID_ST)
            minTree = copy.deepcopy(newTree)
            minCenters_ST = copy.deepcopy(centers_ST)
            minLVCostSum_ST = sumLVCostAtEachStep[len(centers_ST)]
            minTotalCost_ST = newTotalCost_ST

        # Calculate maxDist below for next graph and continue if it is less than 500

        try:  # to check if there is a segment on the graph or there is only one cluster
            minSeg_ST = min(segments_ST, key=lambda obj: obj.getWeight())
            maxDist, tempCenterX, tempCenterY = maxTempInClusterDist(minSeg_ST, clusterByNode_ST, nodesByClusterID_ST)

            if maxDist > distFromT:
                segments_ST.sort(key=lambda obj: obj.getWeight())

                for seg in segments_ST:
                    if seg.getWeight() > distFromT * 2:
                        break
                    else:
                        maxDist, tempCenterX, tempCenterY = maxTempInClusterDist(seg, clusterByNode_ST,
                                                                                 nodesByClusterID_ST)
                        if maxDist <= distFromT:
                            minSeg_ST = seg
                            break
        except:
            break
    outFile.close()
    return minTotalCost_ST, minTree, minCenters_ST, minNodesByClusterID_ST, minLVCostSum_ST


def writeLVDictToText(statsFile, Dict):
    'Writes LVCostDict to a text file for batchPrimsForTransformers.py.'
    outFile = open(statsFile, "w")
    for key in Dict.keys():
        LVCost = Dict[key] * 10
        outFile.write(f"{key} {LVCost} \n")
    outFile.close()
    return 0


def writeCenterSizeToText(statsFile, Dict):
    outFile = open(statsFile, "w")
    for key in Dict.keys():
        size = Dict[key].getWeight()
        outFile.write(f"{size} \n")
    outFile.close()
    return 0


countries = {
  "AF": ["Afghanistan", [
    60.5284298033,
    29.318572496,
    75.1580277851,
    38.4862816432,
  ]],
  "AO": ["Angola", [
    11.6400960629,
    -17.9306364885,
    24.0799052263,
    -4.43802336998,
  ]],
  "AL": ["Albania", [
    19.3044861183,
    39.624997667,
    21.0200403175,
    42.6882473822,
  ]],
  "AE": ["United Arab Emirates", [
    51.5795186705,
    22.4969475367,
    56.3968473651,
    26.055464179,
  ]],
  "AR": ["Argentina", [-73.4154357571, -55.25, -53.628348965, -21.8323104794]],
  "AM": ["Armenia", [
    43.5827458026,
    38.7412014837,
    46.5057198423,
    41.2481285671,
  ]],
  "AQ": ["Antarctica", [-180.0, -90.0, 180.0, -63.2706604895]],
  "TF": ["Fr. S. and Antarctic Lands", [68.72, -49.775, 70.56, -48.625]],
  "AU": ["Australia", [
    113.338953078,
    -43.6345972634,
    153.569469029,
    -10.6681857235,
  ]],
  "AT": ["Austria", [
    9.47996951665,
    46.4318173285,
    16.9796667823,
    49.0390742051,
  ]],
  "AZ": ["Azerbaijan", [
    44.7939896991,
    38.2703775091,
    50.3928210793,
    41.8606751572,
  ]],
  "BI": ["Burundi", [
    29.0249263852,
    -4.49998341229,
    30.752262811,
    -2.34848683025,
  ]],
  "BE": ["Belgium", [
    2.51357303225,
    49.5294835476,
    6.15665815596,
    51.4750237087,
  ]],
  "BJ": ["Benin", [
    0.772335646171,
    6.14215770103,
    3.79711225751,
    12.2356358912,
  ]],
  "BF": ["Burkina Faso", [
    -5.47056494793,
    9.61083486576,
    2.17710778159,
    15.1161577418,
  ]],
  "BD": ["Bangladesh", [
    88.0844222351,
    20.670883287,
    92.6727209818,
    26.4465255803,
  ]],
  "BG": ["Bulgaria", [
    22.3805257504,
    41.2344859889,
    28.5580814959,
    44.2349230007,
  ]],
  "BS": ["Bahamas", [-78.98, 23.71, -77.0, 27.04]],
  "BA": ["Bosnia and Herz.", [15.7500260759, 42.65, 19.59976, 45.2337767604]],
  "BY": ["Belarus", [
    23.1994938494,
    51.3195034857,
    32.6936430193,
    56.1691299506,
  ]],
  "BZ": ["Belize", [
    -89.2291216703,
    15.8869375676,
    -88.1068129138,
    18.4999822047,
  ]],
  "BO": ["Bolivia", [
    -69.5904237535,
    -22.8729187965,
    -57.4983711412,
    -9.76198780685,
  ]],
  "BR": ["Brazil", [
    -73.9872354804,
    -33.7683777809,
    -34.7299934555,
    5.24448639569,
  ]],
  "BN": ["Brunei", [114.204016555, 4.007636827, 115.450710484, 5.44772980389]],
  "BT": ["Bhutan", [
    88.8142484883,
    26.7194029811,
    92.1037117859,
    28.2964385035,
  ]],
  "BW": ["Botswana", [
    19.8954577979,
    -26.8285429827,
    29.4321883481,
    -17.6618156877,
  ]],
  "CF": ["Central African Rep.", [
    14.4594071794,
    2.2676396753,
    27.3742261085,
    11.1423951278,
  ]],
  "CA": ["Canada", [-140.99778, 41.6751050889, -52.6480987209, 83.23324]],
  "CH": ["Switzerland", [
    6.02260949059,
    45.7769477403,
    10.4427014502,
    47.8308275417,
  ]],
  "CL": ["Chile", [-75.6443953112, -55.61183, -66.95992, -17.5800118954]],
  "CN": ["China", [73.6753792663, 18.197700914, 135.026311477, 53.4588044297]],
  "CI": ["Ivory Coast", [
    -8.60288021487,
    4.33828847902,
    -2.56218950033,
    10.5240607772,
  ]],
  "CM": ["Cameroon", [
    8.48881554529,
    1.72767263428,
    16.0128524106,
    12.8593962671,
  ]],
  "CD": ["Congo [Kinshasa]", [
    12.1823368669,
    -13.2572266578,
    31.1741492042,
    5.25608775474,
  ]],
  "CG": ["Congo [Brazzaville]", [
    11.0937728207,
    -5.03798674888,
    18.4530652198,
    3.72819651938,
  ]],
  "CO": ["Colombia", [
    -78.9909352282,
    -4.29818694419,
    -66.8763258531,
    12.4373031682,
  ]],
  "CR": ["Costa Rica", [
    -85.94172543,
    8.22502798099,
    -82.5461962552,
    11.2171192489,
  ]],
  "CU": ["Cuba", [
    -84.9749110583,
    19.8554808619,
    -74.1780248685,
    23.1886107447,
  ]],
  "CY": ["Cyprus", [
    32.2566671079,
    34.5718694118,
    34.0048808123,
    35.1731247015,
  ]],
  "CZ": ["Czech Rep.", [
    12.2401111182,
    48.5553052842,
    18.8531441586,
    51.1172677679,
  ]],
  "DE": ["Germany", [
    5.98865807458,
    47.3024876979,
    15.0169958839,
    54.983104153,
  ]],
  "DJ": ["Djibouti", [41.66176, 10.9268785669, 43.3178524107, 12.6996385767]],
  "DK": ["Denmark", [
    8.08997684086,
    54.8000145534,
    12.6900061378,
    57.730016588,
  ]],
  "DO": ["Dominican Rep.", [
    -71.9451120673,
    17.598564358,
    -68.3179432848,
    19.8849105901,
  ]],
  "DZ": ["Algeria", [
    -8.68439978681,
    19.0573642034,
    11.9995056495,
    37.1183806422,
  ]],
  "EC": ["Ecuador", [
    -80.9677654691,
    -4.95912851321,
    -75.2337227037,
    1.3809237736,
  ]],
  "EG": ["Egypt", [24.70007, 22.0, 36.86623, 31.58568]],
  "ER": ["Eritrea", [36.3231889178, 12.4554157577, 43.0812260272, 17.9983074]],
  "ES": ["Spain", [-9.39288367353, 35.946850084, 3.03948408368, 43.7483377142]],
  "EE": ["Estonia", [
    23.3397953631,
    57.4745283067,
    28.1316992531,
    59.6110903998,
  ]],
  "ET": ["Ethiopia", [32.95418, 3.42206, 47.78942, 14.95943]],
  "FI": ["Finland", [
    20.6455928891,
    59.846373196,
    31.5160921567,
    70.1641930203,
  ]],
  "FJ": ["Fiji", [-180.0, -18.28799, 180.0, -16.0208822567]],
  "FK": ["Falkland Is.", [-61.2, -52.3, -57.75, -51.1]],
  "FR": ["France", [
    -54.5247541978,
    2.05338918702,
    9.56001631027,
    51.1485061713,
  ]],
  "GA": ["Gabon", [
    8.79799563969,
    -3.97882659263,
    14.4254557634,
    2.32675751384,
  ]],
  "GB": ["United Kingdom", [
    -7.57216793459,
    49.959999905,
    1.68153079591,
    58.6350001085,
  ]],
  "GE": ["Georgia", [
    39.9550085793,
    41.0644446885,
    46.6379081561,
    43.553104153,
  ]],
  "GH": ["Ghana", [-3.24437008301, 4.71046214438, 1.0601216976, 11.0983409693]],
  "GN": ["Guinea", [
    -15.1303112452,
    7.3090373804,
    -7.83210038902,
    12.5861829696,
  ]],
  "GM": ["Gambia", [
    -16.8415246241,
    13.1302841252,
    -13.8449633448,
    13.8764918075,
  ]],
  "GW": ["Guinea Bissau", [
    -16.6774519516,
    11.0404116887,
    -13.7004760401,
    12.6281700708,
  ]],
  "GQ": ["Eq. Guinea", [
    9.3056132341,
    1.01011953369,
    11.285078973,
    2.28386607504,
  ]],
  "GR": ["Greece", [
    20.1500159034,
    34.9199876979,
    26.6041955909,
    41.8269046087,
  ]],
  "GL": ["Greenland", [-73.297, 60.03676, -12.20855, 83.64513]],
  "GT": ["Guatemala", [
    -92.2292486234,
    13.7353376327,
    -88.2250227526,
    17.8193260767,
  ]],
  "GY": ["Guyana", [
    -61.4103029039,
    1.26808828369,
    -56.5393857489,
    8.36703481692,
  ]],
  "HN": ["Honduras", [
    -89.3533259753,
    12.9846857772,
    -83.147219001,
    16.0054057886,
  ]],
  "HR": ["Croatia", [13.6569755388, 42.47999136, 19.3904757016, 46.5037509222]],
  "HT": ["Haiti", [
    -74.4580336168,
    18.0309927434,
    -71.6248732164,
    19.9156839055,
  ]],
  "HU": ["Hungary", [
    16.2022982113,
    45.7594811061,
    22.710531447,
    48.6238540716,
  ]],
  "ID": ["Indonesia", [
    95.2930261576,
    -10.3599874813,
    141.03385176,
    5.47982086834,
  ]],
  "IN": ["India", [68.1766451354, 7.96553477623, 97.4025614766, 35.4940095078]],
  "IE": ["Ireland", [
    -9.97708574059,
    51.6693012559,
    -6.03298539878,
    55.1316222195,
  ]],
  "IR": ["Iran", [44.1092252948, 25.0782370061, 63.3166317076, 39.7130026312]],
  "IQ": ["Iraq", [38.7923405291, 29.0990251735, 48.5679712258, 37.3852635768]],
  "IS": ["Iceland", [
    -24.3261840479,
    63.4963829617,
    -13.609732225,
    66.5267923041,
  ]],
  "IL": ["Israel", [
    34.2654333839,
    29.5013261988,
    35.8363969256,
    33.2774264593,
  ]],
  "IT": ["Italy", [6.7499552751, 36.619987291, 18.4802470232, 47.1153931748]],
  "JM": ["Jamaica", [
    -78.3377192858,
    17.7011162379,
    -76.1996585761,
    18.5242184514,
  ]],
  "JO": ["Jordan", [
    34.9226025734,
    29.1974946152,
    39.1954683774,
    33.3786864284,
  ]],
  "JP": ["Japan", [129.408463169, 31.0295791692, 145.543137242, 45.5514834662]],
  "KZ": ["Kazakhstan", [
    46.4664457538,
    40.6623245306,
    87.3599703308,
    55.3852501491,
  ]],
  "KE": ["Kenya", [33.8935689697, -4.67677, 41.8550830926, 5.506]],
  "KG": ["Kyrgyzstan", [
    69.464886916,
    39.2794632025,
    80.2599902689,
    43.2983393418,
  ]],
  "KH": ["Cambodia", [
    102.3480994,
    10.4865436874,
    107.614547968,
    14.5705838078,
  ]],
  "KR": ["S. Korea", [
    126.117397903,
    34.3900458847,
    129.468304478,
    38.6122429469,
  ]],
  "KW": ["Kuwait", [
    46.5687134133,
    28.5260627304,
    48.4160941913,
    30.0590699326,
  ]],
  "LA": ["Laos", [100.115987583, 13.88109101, 107.564525181, 22.4647531194]],
  "LB": ["Lebanon", [
    35.1260526873,
    33.0890400254,
    36.6117501157,
    34.6449140488,
  ]],
  "LR": ["Liberia", [
    -11.4387794662,
    4.35575511313,
    -7.53971513511,
    8.54105520267,
  ]],
  "LY": ["Libya", [9.31941084152, 19.58047, 25.16482, 33.1369957545]],
  "LK": ["Sri Lanka", [
    79.6951668639,
    5.96836985923,
    81.7879590189,
    9.82407766361,
  ]],
  "LS": ["Lesotho", [
    26.9992619158,
    -30.6451058896,
    29.3251664568,
    -28.6475017229,
  ]],
  "LT": ["Lithuania", [
    21.0558004086,
    53.9057022162,
    26.5882792498,
    56.3725283881,
  ]],
  "LU": ["Luxembourg", [
    5.67405195478,
    49.4426671413,
    6.24275109216,
    50.1280516628,
  ]],
  "LV": ["Latvia", [21.0558004086, 55.61510692, 28.1767094256, 57.9701569688]],
  "MA": ["Morocco", [
    -17.0204284327,
    21.4207341578,
    -1.12455115397,
    35.7599881048,
  ]],
  "MD": ["Moldova", [
    26.6193367856,
    45.4882831895,
    30.0246586443,
    48.4671194525,
  ]],
  "MG": ["Madagascar", [
    43.2541870461,
    -25.6014344215,
    50.4765368996,
    -12.0405567359,
  ]],
  "MX": ["Mexico", [-117.12776, 14.5388286402, -86.811982388, 32.72083]],
  "MK": ["Macedonia", [20.46315, 40.8427269557, 22.9523771502, 42.3202595078]],
  "ML": ["Mali", [-12.1707502914, 10.0963607854, 4.27020999514, 24.9745740829]],
  "MM": ["Myanmar", [
    92.3032344909,
    9.93295990645,
    101.180005324,
    28.335945136,
  ]],
  "ME": ["Montenegro", [18.45, 41.87755, 20.3398, 43.52384]],
  "MN": ["Mongolia", [
    87.7512642761,
    41.5974095729,
    119.772823928,
    52.0473660345,
  ]],
  "MZ": ["Mozambique", [
    30.1794812355,
    -26.7421916643,
    40.7754752948,
    -10.3170960425,
  ]],
  "MR": ["Mauritania", [
    -17.0634232243,
    14.6168342147,
    -4.92333736817,
    27.3957441269,
  ]],
  "MW": ["Malawi", [
    32.6881653175,
    -16.8012997372,
    35.7719047381,
    -9.23059905359,
  ]],
  "MY": ["Malaysia", [
    100.085756871,
    0.773131415201,
    119.181903925,
    6.92805288332,
  ]],
  "NA": ["Namibia", [
    11.7341988461,
    -29.045461928,
    25.0844433937,
    -16.9413428687,
  ]],
  "NC": ["New Caledonia", [
    164.029605748,
    -22.3999760881,
    167.120011428,
    -20.1056458473,
  ]],
  "NE": ["Niger", [
    0.295646396495,
    11.6601671412,
    15.9032466977,
    23.4716684026,
  ]],
  "NG": ["Nigeria", [
    2.69170169436,
    4.24059418377,
    14.5771777686,
    13.8659239771,
  ]],
  "NI": ["Nicaragua", [
    -87.6684934151,
    10.7268390975,
    -83.147219001,
    15.0162671981,
  ]],
  "NL": ["Netherlands", [
    3.31497114423,
    50.803721015,
    7.09205325687,
    53.5104033474,
  ]],
  "NO": ["Norway", [4.99207807783, 58.0788841824, 31.29341841, 80.6571442736]],
  "NP": ["Nepal", [80.0884245137, 26.3978980576, 88.1748043151, 30.4227169866]],
  "NZ": ["New Zealand", [
    166.509144322,
    -46.641235447,
    178.517093541,
    -34.4506617165,
  ]],
  "OM": ["Oman", [52.0000098, 16.6510511337, 59.8080603372, 26.3959343531]],
  "PK": ["Pakistan", [
    60.8742484882,
    23.6919650335,
    77.8374507995,
    37.1330309108,
  ]],
  "PA": ["Panama", [
    -82.9657830472,
    7.2205414901,
    -77.2425664944,
    9.61161001224,
  ]],
  "PE": ["Peru", [
    -81.4109425524,
    -18.3479753557,
    -68.6650797187,
    -0.0572054988649,
  ]],
  "PH": ["Philippines", [
    117.17427453,
    5.58100332277,
    126.537423944,
    18.5052273625,
  ]],
  "PG": ["Papua New Guinea", [
    141.000210403,
    -10.6524760881,
    156.019965448,
    -2.50000212973,
  ]],
  "PL": ["Poland", [
    14.0745211117,
    49.0273953314,
    24.0299857927,
    54.8515359564,
  ]],
  "PR": ["Puerto Rico", [
    -67.2424275377,
    17.946553453,
    -65.5910037909,
    18.5206011011,
  ]],
  "KP": ["N. Korea", [
    124.265624628,
    37.669070543,
    130.780007359,
    42.9853868678,
  ]],
  "PT": ["Portugal", [
    -9.52657060387,
    36.838268541,
    -6.3890876937,
    42.280468655,
  ]],
  "PY": ["Paraguay", [
    -62.6850571357,
    -27.5484990374,
    -54.2929595608,
    -19.3427466773,
  ]],
  "QA": ["Qatar", [50.7439107603, 24.5563308782, 51.6067004738, 26.1145820175]],
  "RO": ["Romania", [20.2201924985, 43.6884447292, 29.62654341, 48.2208812526]],
  "RU": ["Russia", [-180.0, 41.151416124, 180.0, 81.2504]],
  "RW": ["Rwanda", [
    29.0249263852,
    -2.91785776125,
    30.8161348813,
    -1.13465911215,
  ]],
  "SA": ["Saudi Arabia", [
    34.6323360532,
    16.3478913436,
    55.6666593769,
    32.161008816,
  ]],
  "SD": ["Sudan", [21.93681, 8.61972971293, 38.4100899595, 22.0]],
  "SS": ["S. Sudan", [23.8869795809, 3.50917, 35.2980071182, 12.2480077571]],
  "SN": ["Senegal", [
    -17.6250426905,
    12.332089952,
    -11.4678991358,
    16.5982636581,
  ]],
  "SB": ["Solomon Is.", [
    156.491357864,
    -10.8263672828,
    162.398645868,
    -6.59933847415,
  ]],
  "SL": ["Sierra Leone", [
    -13.2465502588,
    6.78591685631,
    -10.2300935531,
    10.0469839543,
  ]],
  "SV": ["El Salvador", [
    -90.0955545723,
    13.1490168319,
    -87.7235029772,
    14.4241327987,
  ]],
  "SO": ["Somalia", [40.98105, -1.68325, 51.13387, 12.02464]],
  "RS": ["Serbia", [18.82982, 42.2452243971, 22.9860185076, 46.1717298447]],
  "SR": ["Suriname", [
    -58.0446943834,
    1.81766714112,
    -53.9580446031,
    6.0252914494,
  ]],
  "SK": ["Slovakia", [
    16.8799829444,
    47.7584288601,
    22.5581376482,
    49.5715740017,
  ]],
  "SI": ["Slovenia", [
    13.6981099789,
    45.4523163926,
    16.5648083839,
    46.8523859727,
  ]],
  "SE": ["Sweden", [
    11.0273686052,
    55.3617373725,
    23.9033785336,
    69.1062472602,
  ]],
  "SZ": ["Swaziland", [
    30.6766085141,
    -27.2858794085,
    32.0716654803,
    -25.660190525,
  ]],
  "SY": ["Syria", [35.7007979673, 32.312937527, 42.3495910988, 37.2298725449]],
  "TD": ["Chad", [13.5403935076, 7.42192454674, 23.88689, 23.40972]],
  "TG": ["Togo", [
    -0.0497847151599,
    5.92883738853,
    1.86524051271,
    11.0186817489,
  ]],
  "TH": ["Thailand", [
    97.3758964376,
    5.69138418215,
    105.589038527,
    20.4178496363,
  ]],
  "TJ": ["Tajikistan", [
    67.4422196796,
    36.7381712916,
    74.9800024759,
    40.9602133245,
  ]],
  "TM": ["Turkmenistan", [
    52.5024597512,
    35.2706639674,
    66.5461503437,
    42.7515510117,
  ]],
  "TL": ["East Timor", [
    124.968682489,
    -9.39317310958,
    127.335928176,
    -8.27334482181,
  ]],
  "TT": ["Trinidad and Tobago", [-61.95, 10.0, -60.895, 10.89]],
  "TN": ["Tunisia", [
    7.52448164229,
    30.3075560572,
    11.4887874691,
    37.3499944118,
  ]],
  "TR": ["Turkey", [
    26.0433512713,
    35.8215347357,
    44.7939896991,
    42.1414848903,
  ]],
  "TW": ["Taiwan", [
    120.106188593,
    21.9705713974,
    121.951243931,
    25.2954588893,
  ]],
  "TZ": ["Tanzania", [29.3399975929, -11.7209380022, 40.31659, -0.95]],
  "UG": ["Uganda", [29.5794661801, -1.44332244223, 35.03599, 4.24988494736]],
  "UA": ["Ukraine", [
    22.0856083513,
    44.3614785833,
    40.0807890155,
    52.3350745713,
  ]],
  "UY": ["Uruguay", [
    -58.4270741441,
    -34.9526465797,
    -53.209588996,
    -30.1096863746,
  ]],
  "US": ["United States", [-171.791110603, 18.91619, -66.96466, 71.3577635769]],
  "UZ": ["Uzbekistan", [
    55.9289172707,
    37.1449940049,
    73.055417108,
    45.5868043076,
  ]],
  "VE": ["Venezuela", [
    -73.3049515449,
    0.724452215982,
    -59.7582848782,
    12.1623070337,
  ]],
  "VN": ["Vietnam", [
    102.170435826,
    8.59975962975,
    109.33526981,
    23.3520633001,
  ]],
  "VU": ["Vanuatu", [
    166.629136998,
    -16.5978496233,
    167.844876744,
    -14.6264970842,
  ]],
  "PS": ["West Bank", [
    34.9274084816,
    31.3534353704,
    35.5456653175,
    32.5325106878,
  ]],
  "YE": ["Yemen", [42.6048726743, 12.5859504257, 53.1085726255, 19.0000033635]],
  "ZA": ["South Africa", [
    16.3449768409,
    -34.8191663551,
    32.830120477,
    -22.0913127581,
  ]],
  "ZM": ["Zambia", [
    21.887842645,
    -17.9612289364,
    33.4856876971,
    -8.23825652429,
  ]],
  "ZW": ["Zimbabwe", [
    25.2642257016,
    -22.2716118303,
    32.8498608742,
    -15.5077869605,
  ]],
};

def generateDictsFromShp(x, y, pue):

    nodesByClusterID = collections.defaultdict(list)
    clusterByNode = {}
    nodes = {}
    centers = {}
    LVCostDict = {}
    nodes_weights_output = []
    FID = 0
    for xx, yy, pp in zip(x, y, pue):
        if pp:
            nodeWeight = 10
        else:
            nodeWeight = 1
        FID += 1
        nodes[FID] = Node(FID, xx, yy, nodeWeight)
        centers[FID] = Node(FID, xx, yy, nodeWeight)
        clusterByNode[nodes[FID]] = FID
        nodesByClusterID[FID].append(nodes[FID])
        LVCostDict[FID] = 0
        nodes_weights_output.append(nodeWeight)

    return nodesByClusterID, clusterByNode, nodes, centers, LVCostDict ,nodes_weights_output



def tlnd(outputDir, df, pues_in_cell=0, structures_in_cell=0):
    tf = os.path.join(outputDir, "transformers.csv")
    if os.path.exists(tf):
        with open(tf) as tf_file:
            raw = list(csv.DictReader(tf_file))
            
            latitudes = [d['latitude'] for d in raw]
            longitudes = [d['longitude'] for d in raw]
            lv_per_customer = [d['lv_meters_per_customer'] for d in raw]
            customers = [d['customers'] for d in raw]
            pues = [d['pues'] for d in raw]

            output = {
                "latitude": latitudes,
                "longitude": longitudes,
                "lv_meters_per_customer": lv_per_customer,
                "customers": customers,
                "pues": pues
            }
            return output
        
    searchRadius = 100000  # meters
    # Cost parameters:
    MV = 25  # Cost of MV per meter
    LV = 10  # Cost of LV per meter
    TCost = 2000  # Transformer Cost
    distFromT = 750  # Dmax, direct distance from transformers
    maxLVLenghtInCluster = 1000  # Lmax
    logfilename = outputDir + '/' + 'modelStatus.txt'
    startTime = time.time()

    zone_number = df.utm_zone_number[0]
    zone_letter = df.utm_zone_letter[0]
    #print("UTM", zone_number, zone_letter)
    x = df.x
    y = df.y
    pue = df.pue
    nodesByClusterID, clusterByNode, nodes, centers, LVCostDict , _ = generateDictsFromShp(x, y, pue)

    _, tree, centers, nodesByClusterID, _ = run(centers, nodesByClusterID, clusterByNode,
                                                                    LVCostDict, searchRadius, MV, LV, TCost,
                                                                    distFromT,
                                                                    maxLVLenghtInCluster, outputDir, logfilename)

    statsFile1 = outputDir + os.sep + "LVCostDict.txt"
    statsFile2 = outputDir + os.sep + "CenterSize.txt"
    writeLVDictToText(statsFile1, LVCostDict)
    writeCenterSizeToText(statsFile2, centers)
    MVLength = tree.getTotalEdgeWeight()
    MVCost = MVLength * MV
    numTransformer = len(centers)

    try:
        netID = tree.getNetID(centers.values()[0])
    except:
        netID = 0
        tree._nodesByNetID[0] = []
        tree._network[netID] = []
    
    my_lv = 0
    transformers = []

    for ID in centers.keys():
        nodesByNodeID = {}
        # Start on -1 to not over count transformers
        customers = -1
        totalWeight = 0
        segments, lvCost = CMST(nodesByClusterID[ID], maxLVLenghtInCluster, centers[ID])

        my_lv += lvCost
        for segment in segments.values():
            node1 = segment.getNode1()
            node2 = segment.getNode2()
            if node1.getID() not in nodesByNodeID.keys():
                nodesByNodeID[node1.getID()] = node1
            if node2.getID() not in nodesByNodeID.keys():
                nodesByNodeID[node2.getID()] = node2

        for node in nodesByNodeID.values():
            tree._netIDByNode[node] = netID
            tree._nodesByNetID[netID].append(node)
            customers += 1
            totalWeight += node._weight
        for segment in segments.values():
            tree._network[netID].append(segment)

        pues = int((totalWeight - customers)/10)
        transformers.append((centers[ID]._x, centers[ID]._y, lvCost/customers, customers, pues))
   
    
    lat, lon = to_latlon(centers[ID]._x, centers[ID]._y, zone_number, zone_letter)

    coords = [to_latlon(dd[0], dd[1], zone_number, zone_letter) for dd in transformers]

    #coords = [(dd[0], dd[1]) for dd in transformers]
    latitudes = [cc[1] for cc in coords]
    longitudes = [cc[0] for cc in coords]

    lv_per_customer = [dd[2] for dd in transformers]
    customers = [dd[3] for dd in transformers]
    pues = [dd[4] for dd in transformers]
    
    output = {"latitude": latitudes, "longitude": longitudes,
              "lv_meters_per_customer": lv_per_customer,
              "customers": customers, "pues": pues}
    
    rows = list(zip(output["latitude"], output["longitude"], 
                output["lv_meters_per_customer"], output["customers"], output["pues"]))

    with open(os.path.join(outputDir, "transformers.csv"), "w", newline="") as f:
        csvwriter = csv.writer(f)

        # Write headers
        csvwriter.writerow(output.keys())

        # Write data
        csvwriter.writerows(rows)

    with open(outputDir + "/" + 'modelOutput.txt', 'w') as dst:
        dst.write("NumStructures:" + str(len(nodes)) + "\n")
        dst.write("LVLength:" + str(my_lv) + "\n")
        dst.write("LVPerCustomer:" + str(float(my_lv) / len(nodes)) + "\n")
        dst.write("MVLength:" + str(MVLength) + "\n")
        dst.write("MVPerCustomer:" + str(MVLength / len(nodes)) + "\n")
        dst.write("Num Transformers:" + str(numTransformer) + "\n")
        dst.write("Customers Per Tx:" + str(len(nodes) / float(numTransformer)) + "\n")
        dst.write("Total LV Cost:" + str(my_lv * float(LV)) + "\n")
        dst.write("Total MV Cost:" + str(MVCost) + "\n")
        dst.write("StructuresInCell:" + str(structures_in_cell) + "\n")
        dst.write("PuesInCell:" + str(pues_in_cell) + "\n")
        transformerCost = numTransformer * TCost
        dst.write("Transformer Cost:" + str(transformerCost) + "\n")
        total_cost = MVCost + my_lv * float(LV) + transformerCost
        dst.write("Total Cost:" + str(total_cost) + "\n")
        runningT = time.time() - startTime
        dst.write("Total Running Time:" + str(runningT) + "\n")
        # with open(outputDir+'modelOutput.txt', 'a') as dst:
        runningT1 = time.time() - startTime
        dst.write("Final Running Time:" + str(runningT1))


    return output
  


def read_lan(path_to_gzipped_csv):
    with gzip.open(path_to_gzipped_csv, 'rt') as f:
        data = list(csv.DictReader(f))
    
    x = [float(dd["x"]) for dd in data]
    y = [float(dd["y"]) for dd in data]
    
    return x, y

def start():
    import geopandas as gpd
    import pandas as pd
    import osgeo.osr as osr
    import duckdb as db

    srs = osr.SpatialReference()
    srs.SetFromUserInput("EPSG:32636")
    uuid = sys.argv[2]
    input_file = sys.argv[3]
    out = sys.argv[4]
    output_dir = os.path.join(out, uuid)

    structures_raw = db.sql(f"SELECT lan, structure, pue FROM '{input_file}/*.parquet'")
    structures = db.sql(f"SELECT structure, pue from structures_raw WHERE lan = {uuid}").df()
    h3_structures = [str(hex(ss))[2:] for ss in structures.structure]

    geo = [h3.h3_to_geo(h3_s) for h3_s in h3_structures] 
    x_y = [from_latlon(latitude=ll[1], longitude=ll[0]) for ll in geo]
   
    xy_df = pd.DataFrame(x_y, columns=['x', 'y', 'utm_zone_number', 'utm_zone_letter'])
    
    df = pd.concat([structures, xy_df], axis=1)

    non_pue = df.pue.isna().sum()
    structures_in_cell = len(df)
    pues_in_cell = structures_in_cell - non_pue

    print(uuid, " pues ", pues_in_cell)
    print(uuid, " structures ", structures_in_cell)

    if len(df) < structures_in_cell:
        assert False, "How this did happen?"

    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    
    if len(df) > 0:
        tlnd(output_dir, df, pues_in_cell=pues_in_cell, structures_in_cell=structures_in_cell)


if __name__ == "__main__":
    start()
    #import cProfile
    #cProfile.run('start()')

