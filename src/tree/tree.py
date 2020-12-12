#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 17:19:12 2020

@author: leonard
"""
from Constants import BITS_FOR_POSITIONS, MAX_INT, TREE_NUM_BEFORE_NODESPLIT
import numpy as np
from sys import exit

class ngbnode():
    "A structure for node data"
    def __init__(self):
        #By default initialize the root node
        self.center   =  2**(BITS_FOR_POSITIONS - 1)   #The center of the node
        self.level    =  0                             #A useful way to store the sidelength of the node
        
        #refernce information to other nodes
        self.sibling  = -1                          
        self.nextnode = -1
        self.father   = -1
        #Is this node empty?
        self.notEmpty = 0
        
        #SPH data
        self.Mass              = 0.0
        self.center_offset_min = 0
        self.center_offset_max = 0
        
class ngbtree():
    "A neighbor search tree"
    def __init__(self, particles):
        #List of particles
        N = len(particles)
        if N == 0:
            print("There are no particles!\n")
            exit()
        
        self.Tp          = particles                     

        self.Father      = np.empty(N, dtype = int)
        self.Nextnode    = np.empty(N + 2, dtype = int)
        self.Nodes       = list()
        #np.empty(N + 1 + N + 100, dtype = ngbnode)
        
        #data for referencing nodes/particles
        self.MaxPart              = N
        self.MaxNodes             = 1 + N + 100
        self.NextFreeNode         = 0
        self.FirstNonTopLevelNode = 0
        
        #Now build the tree
        self.treebuild()
        
    def get_nodep(self, no):
        return self.Nodes[no-self.MaxPart]
        
    def treebuild(self):
        "construct the tree and fill the nodes with all the necessary information"
        self.treeconstruct()
        #start updating the tree nodes starting from the root node
        self.update_node_recursive(self.MaxPart)
            
        
    def update_node_recursive(self, no):
        "walk down the tree to update the SPH information in all nodes"
        nop = self.get_nodep(no)
        #initialise the SPH info
        mass      = 0.0
        not_empty = 0
        halflen = 2**(BITS_FOR_POSITIONS - 1 - nop.level)
        center_offset_min = halflen -1
        center_offset_max = -halflen
            
        p = nop.nextnode
        while p != nop.sibling:
            if p>=0:
                if p >= self.MaxPart:
                    #node
                    self.update_node_recursive(p)
                if p < self.MaxPart:
                    #particle
                    mass += self.Tp[p].mass
                        
                    offset = self.Tp[p].position - nop.center
                    if offset < center_offset_min:
                        center_offset_min = offset
                    if offset > center_offset_max:
                        center_offset_max = offset
                        
                    not_empty = 1
                    p = self.Nextnode[p]
                else:
                    #node
                    node = self.get_nodep(p)
                    
                    mass += node.Mass
                        
                    offset_min = node.center_offset_min + node.center \
                                                        - nop.center
                    offset_max = node.center_offset_max + node.center \
                                                        - nop.center
                                                            
                    if offset_min < center_offset_min:
                        center_offset_min = offset_min
                    if offset_max > center_offset_max:
                        center_offset_max = offset_max
                        
                    not_empty |= node.notEmpty
                    p = node.sibling
        #now update the node
        self.get_nodep(no).Mass              = mass
        self.get_nodep(no).notEmpty          = not_empty
        self.get_nodep(no).center_offset_min = center_offset_min
        self.get_nodep(no).center_offset_max = center_offset_max 
        
    def treeconstruct(self):
        "create the root node and recursively all the daughternodes and so on"
        #create the root node
        self.NextFreeNode = self.MaxPart
        self.Nodes.append(ngbnode())
        self.NextFreeNode +=1
            
        index_list = list()
        for i in range(self.MaxPart):
            index_list.append([i, 0])
        self.insert_points(self.MaxPart, index_list, self.MaxPart, 0, -1)
        
    def insert_points(self, num, index_list, th, level, sibling):
        "Recursively constructs the tree nodes and fills them with particles"
        if level >= BITS_FOR_POSITIONS:
            print("Need higher resolution! BITS_FOR_POSITIONS = %d"\
                  %(BITS_FOR_POSITIONS))
            exit()
        mask = 1 << (BITS_FOR_POSITIONS - 1 - level)
        shift = (BITS_FOR_POSITIONS - 1 - level)
        centermask = ~(MAX_INT >> level) & MAX_INT
        subcount  = [0, 0]
        subnode   = [0, 0]
        subintpos = [0, 0]
        #determine the nodes in which the particles will be put
        for i in range(num):
            p = index_list[i][0]
            intpos = self.Tp[p].position
            snode = (intpos & mask) >> shift
            subcount[snode] += 1
            subintpos[snode] = intpos
            index_list[i][1] = snode
                
        centermask >>= 1
        centermask |= ~(MAX_INT >> 1) & MAX_INT
        mask >>= 1
            
        #create the daughter nodes
        nfreep_last = 0
        for i in range(2):
            if subcount[i] > TREE_NUM_BEFORE_NODESPLIT:
                thnew = self.NextFreeNode
                self.NextFreeNode +=1
                subnode[i] = thnew
                if nfreep_last:
                    self.get_nodep(thnew - 1).sibling  = thnew
                    self.get_nodep(thnew - 1).nextnode = thnew
                else:
                    self.get_nodep(th).nextnode = thnew
                nfreep = ngbnode()
                nfreep.father = th
                nfreep.level  = level + 1
                nfreep.center = (subintpos[i] & centermask) | mask
                self.Nodes.append(nfreep)
                nfreep_last = nfreep
            
        # now insert the particles that are chained to the node
        p_last = -1
        p_first = -1
        off = 0
        # now insert the particle groups into the created daughter nodes, or chain them to the node
        for count in subcount:
            if count <= TREE_NUM_BEFORE_NODESPLIT:
                # put the particles into the node as a chain
                for i in range(count):
                    p = index_list[off + i][0]
                    if not nfreep_last and p_first == -1:
                        self.get_nodep(th).nextnode = p
                        
                    self.Father[p] = th
                       
                    if p_last >=0:
                        self.Nextnode[p_last] = p
                    p_last = p
                    if p_first == -1:
                        p_first = p
            off += count   
            
        if p_last >= 0:
            self.Nextnode[p_last] = sibling
            
        if nfreep_last:
            if p_last >= 0:
                self.get_nodep(thnew).sibling  = p_first
                self.get_nodep(thnew).nextnode = p_first
            else:
                self.get_nodep(thnew).sibling  = sibling
                self.get_nodep(thnew).nextnode = sibling
            
        off = 0
        for i in range(2):
            if subcount[i] > TREE_NUM_BEFORE_NODESPLIT:
                if subnode[i] < self.MaxPart + 1:
                    print("subnode[i]=%d < MaxPart=%d + 1" \
                              %(subnode[i], self.MaxPart))
                    exit()
                self.insert_points(subcount[i], \
                                   index_list[off : off + subcount[i]], \
                                   subnode[i], level + 1, \
                                   self.get_nodep(subnode[i]).sibling)
            off += subcount[i]