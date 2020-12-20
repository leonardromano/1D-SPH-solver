#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 11:19:30 2020

@author: leonard
"""
from Constants import K, DK, Mp, FacIntToCoord, BITS_FOR_POSITIONS, DIM, \
                      BoundaryPeriodic, NORM_COEFF
from src.sph.Kernel import kernel 
from src.data.Particle_Class import Particle
from src.data.int_conversion import get_distance_vector
import numpy as np
from math import ceil
from sys import exit
    
def sph_density_evaluate_particle_node_opening_criterion(particle, node):
    """
    This function checks whether there is a spatial overlap between the 
    (rectangular) enclosing box of the particles contained in a node, 
    and the search region.
    """
    if node.level <= 0:
        return 1
    left  = get_distance_vector(node.center_offset_min + node.center, \
                                particle.position \
                                - ceil(particle.smoothingLength/FacIntToCoord))
    right = get_distance_vector(node.center_offset_max + node.center, \
                                particle.position \
                                - ceil(particle.smoothingLength/FacIntToCoord))
    for i in range(DIM):
        if left[i] > 2 * ceil(particle.smoothingLength/FacIntToCoord) and \
           right[i] > left[i]:
               return 0
    return 1

def sph_density_open_node(particle, nop, NgbTree, ahead = False):
    "Continues to walk the tree for the particle by opening a node."
    p = nop.nextnode
    while p != nop.sibling:
        if p < 0:
            print("p=%d < 0  node.sibling=%d node.nextnode=%d" \
                   %(p, nop.sibling, nop.nextnode))
            exit()
        nextp = 0
        typep = ""
        if p < NgbTree.MaxPart:
            nextp = NgbTree.Nextnode[p]
            typep = "particle"
        else:
            node = NgbTree.get_nodep(p)
            nextp = node.sibling
            typep = "node"
        sph_density_interact(particle, p, typep, NgbTree, ahead)
        
        p = nextp

def sph_density_interact(particle, no, no_type, NgbTree, ahead = False):
    """
    Take care of SPH density interaction between the particle, and the node
    referenced through no. The node can either be a node or a particle.
    """
    if no_type == "particle":
        #we have a particle check whether it's a neighbor
        ngb = NgbTree.Tp[no]
        dx = get_distance_vector(particle.position, ngb.position)
        if np.linalg.norm(dx) > ceil(particle.smoothingLength/FacIntToCoord):
            return
        #compute relative velocity
        if ahead:
            vij = particle.velocity_ahead - ngb.velocity_ahead 
        else:
            vij = particle.velocity - ngb.velocity
        
        particle.neighbors.append([ngb, dx, vij])
    else:
        node = NgbTree.get_nodep(no)
        if not node.notEmpty:
            return
        openflag = sph_density_evaluate_particle_node_opening_criterion(particle, node)
        if openflag:
            sph_density_open_node(particle, node, NgbTree, ahead)

def densities_determine(targetlist, NgbTree, particles, ahead = False):
    """
    for each target walk the tree to determine the neighbors and then 
    compute density and thermodynamic quantities
    """
    for i in targetlist:
        particle = particles[i]
        particle.density            = 0
        particle.dhsmlDensityFactor = 0
        particle.neighbors = list()
        sph_density_interact(particle, NgbTree.MaxPart, "node", NgbTree, ahead)
        # If the particle is close to the wall we need to add ghosts
        add_ghosts(particle, ceil(particle.smoothingLength/FacIntToCoord), ahead)
        evaluate_kernel(particle)
        set_thermodynamic_variables(particle, ahead)

def density(particles, NgbTree, ahead = False):
    "For each particle compute density, smoothing length and thermodynamic variables"
    Left = np.zeros(len(particles))
    Right = np.zeros(len(particles))
    targetlist = list()
    for particle in particles:
        targetlist.append(particle.index)
    while True:
        # now do the primary work with this call
        densities_determine(targetlist, NgbTree, particles, ahead)
        # do final operations on results
        npleft = 0
        for p in range(len(targetlist)):
            i = targetlist[p]
            particle = particles[i]
            numberOfNeighbors = NORM_COEFF * particle.smoothingLength**(DIM) \
                                * particle.density/Mp
            if (abs(numberOfNeighbors-K)>DK):    
                #check whether we're done
                if Left[i]>0 and Right[i] >0:
                    if Right[i]-Left[i] < 10e-3 * Left[i]:
                        continue
                #need to redo this particle
                targetlist[npleft] = i
                npleft += 1
                Left[i], Right[i] = update_bounds(Left[i], Right[i], \
                                                  numberOfNeighbors, \
                                                  particle.smoothingLength)
                update_smoothing_length(Left[i], Right[i], \
                                        numberOfNeighbors, particle)
        targetlist = targetlist[:npleft]
        if len(targetlist) <=0:
            break
    # Need to initialize the ghost particles thermodynamic quantities
    for particle in particles:
        update_ghosts(particle, particles, ceil(particle.smoothingLength/FacIntToCoord), ahead)
    
def update_bounds(lowerBound, upperBound, numberOfNeighbors, h):
    "Update the bounds for the smoothing length in the bisection algorithm"
    if numberOfNeighbors < K-DK:
        lowerBound = max(lowerBound, h)
    else:
        if upperBound != 0:
            if h < upperBound:
                upperBound = h
        else:
            upperBound = h
    return lowerBound, upperBound

def update_smoothing_length(lowerBound, upperBound, numberOfNeighbors, particle):
    "Perform the bisection part of the bisection algorithm"
    if lowerBound >0 and upperBound >0:
        particle.smoothingLength = ((lowerBound**3 + upperBound**3)/2)**(1/3)
    else:
        if upperBound == 0 and lowerBound == 0:
            print("Upper and Lower bounds not updated!")
            exit()
            
        if upperBound == 0 and lowerBound>0:
            if abs(numberOfNeighbors - K) < 0.5*K:
                fac = 1-(numberOfNeighbors - K)/numberOfNeighbors/DIM*\
                    particle.dhsmlDensityFactor
                if fac<1.26:
                    particle.smoothingLength *=fac
                else:
                    particle.smoothingLength *=1.26
            else:
                particle.smoothingLength *=1.26
                    
        if upperBound > 0 and lowerBound  == 0:
            if abs(numberOfNeighbors - K) < 0.5*K:
                fac = 1-(numberOfNeighbors - K)/numberOfNeighbors/DIM*\
                    particle.dhsmlDensityFactor
                if fac > 1/1.26:
                    particle.smoothingLength *=fac
                else:
                    particle.smoothingLength /=1.26
            else:
                particle.smoothingLength /=1.26
                

def add_ghosts(particle, intHsml, ahead = False):
    for i in range(DIM):
        if not BoundaryPeriodic[i]:
            dist_from_wall = get_distance_from_wall(particle, i)
            if dist_from_wall < intHsml:
                #this particle is close to the wall
                min_dist_from_wall = \
                    get_minimum_distance_from_wall(particle.neighbors, i)
                listOfNeighbors = copy_list(particle.neighbors)
                for [neighbor, intDistance, vij] in listOfNeighbors:
                    ngb_dist_from_wall = get_distance_from_wall(neighbor, i)
                    if ngb_dist_from_wall > dist_from_wall and \
                       dist_from_wall < abs(intDistance[i]) + \
                       0.75 * min_dist_from_wall:
                           particle.neighbors.\
                               append(mirror_particle(particle, neighbor, \
                                                      i, ahead))

def update_ghosts(particle, particles, intHsml, ahead = False):
    min_dist_from_wall = (1 << BITS_FOR_POSITIONS)
    for i in range(DIM):
        dist_from_wall = get_distance_from_wall(particle, i)
        if dist_from_wall < min_dist_from_wall:
            min_dist_from_wall = dist_from_wall
    if min_dist_from_wall < intHsml:
        for [neighbor, intDistance, vij] in particle.neighbors:
            if neighbor.density == 0:
                #we have a mirror particle with empty SPH data
                neighbor.density = particles[neighbor.index].density
                neighbor.dhsmlDensityFactor = particles[neighbor.index].dhsmlDensityFactor
                neighbor.update_pressure(ahead)
                neighbor.update_soundspeed()
                
    
def evaluate_kernel(particle):
    "Perform the neighbor sum to compute density and SPH correction factor"
    for [neighbor, intDistance, vij] in particle.neighbors:
        #precompute the kernel and its derivative
        dist = np.linalg.norm(intDistance) * FacIntToCoord
        kern = kernel(dist/particle.smoothingLength, particle.smoothingLength)
        dkern = kernel(dist/particle.smoothingLength, particle.smoothingLength, \
                       True)
        #add neighbor's contribution to density and dhsml
        particle.density += neighbor.mass * kern
        particle.dhsmlDensityFactor -= neighbor.mass * dist * dkern
    
    #now if we have more than one neighbor update the dhsml factor
    if particle.dhsmlDensityFactor > 0:
        particle.dhsmlDensityFactor = DIM * particle.density/particle.dhsmlDensityFactor
        if particle.dhsmlDensityFactor > 10:
            print("Number of neighbors = %d, index = %d"\
              %(len(particle.neighbors), particle.index))
    else:
        print("Number of neighbors = %d, index = %d"\
              %(len(particle.neighbors), particle.index))
        particle.dhsmlDensityFactor = 1
        
def set_thermodynamic_variables(particle, ahead = False):
    "Update the particle's pressure and sound speed"
    particle.update_pressure(ahead)
    particle.update_soundspeed()

def get_distance_from_wall(particle, axis):
    "returns the distance from the nearest non-periodic boundary"
    return int(min(particle.position[axis], \
                   (1 << BITS_FOR_POSITIONS) - particle.position[axis]))

def get_minimum_distance_from_wall(particles, axis):
    "Determines the minimum distance from the wall among a list of particles"
    minimumDistanceFromWall = 1 << BITS_FOR_POSITIONS
    for particle in particles:
        dist = get_distance_from_wall(particle[0], axis)
        if dist < minimumDistanceFromWall:
            minimumDistanceFromWall = dist
    return minimumDistanceFromWall
        
def mirror_particle(particle, neighbor, axis, ahead):
    """
    creates an instance of a ghost particle mirroring a neighbor 
    of a particle close to the wall
    """
    velocity = np.zeros(DIM, dtype = float)
    if ahead: 
        velocity += 2 * particle.velocity_ahead - neighbor.velocity_ahead
        vij = particle.velocity_ahead - velocity
    else:
        velocity += 2 * particle.velocity - neighbor.velocity
        vij = particle.velocity - velocity
    position  = np.zeros(DIM, dtype = int)
    position += neighbor.position
    position[axis] += 2*(particle.position[axis] - neighbor.position[axis])
    mirrorParticle = Particle(position, \
                              velocity, neighbor.entropy, neighbor.index)
    mirrorParticle.entropy = neighbor.entropy
    mirrorParticle.entropy_ahead = neighbor.entropy_ahead
    mirrorParticle.smoothingLength = neighbor.smoothingLength
    return [mirrorParticle, particle.position - mirrorParticle.position, vij]

def copy_list(oldList):
    "creates a hard copy of a list"
    newList = list()
    for item in oldList:
        newList.append(item)
    return newList