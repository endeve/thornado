#!/usr/bin/env python

import sys
import os

import glob
import re

import numpy as np

import h5py as h5

class thornado_hdf5:

    filename = {}
    path     = {}

    time_dsets  = {
                   't':'/Time'
                  }

    grid_dsets  = {
                   'X1'  :'/Spatial Grid/X1',
                   'X2'  :'/Spatial Grid/X2',
                   'X3'  :'/Spatial Grid/X3',
                   'X1_C':'/Spatial Grid/X1_C',
                   'X2_C':'/Spatial Grid/X2_C',
                   'X3_C':'/Spatial Grid/X3_C',
                   'dX1' :'/Spatial Grid/dX1',
                   'dX2' :'/Spatial Grid/dX2',
                   'dX3' :'/Spatial Grid/dX3'
                  }

    energy_dsets = {
                    'E'  :'/Energy Grid/E',
                    'E_C':'/Energy Grid/E_C',
                    'dE' :'/Energy Grid/dE'
                   }

    fluid_dsets = {
                   'dens'  :'/Fluid Fields/Conserved/Conserved Baryon Density',
                   'densYe':'/Fluid Fields/Conserved/Conserved Electron Density',
                   'tau'   :'/Fluid Fields/Conserved/Conserved Energy Density',
                   'scon1' :'/Fluid Fields/Conserved/Conserved Momentum Density (1)',
                   'scon2' :'/Fluid Fields/Conserved/Conserved Momentum Density (2)',
                   'scon3' :'/Fluid Fields/Conserved/Conserved Momentum Density (3)',
                   'rho'   :'/Fluid Fields/Primitive/Comoving Baryon Density',
                   'rhoYe' :'/Fluid Fields/Primitive/Comoving Electron Density',
                   'eps'   :'/Fluid Fields/Primitive/Internal Energy Density',
                   'vel1'  :'/Fluid Fields/Primitive/Three-Velocity (1)',
                   'vel2'  :'/Fluid Fields/Primitive/Three-Velocity (2)',
                   'vel3'  :'/Fluid Fields/Primitive/Three-Velocity (3)',
                   'xAlpha':'/Fluid Fields/Auxiliary/Alpha Mass Fraction',
                   'mu_e'  :'/Fluid Fields/Auxiliary/Electron Chemical Potential',
                   'Ye'    :'/Fluid Fields/Auxiliary/Electron Fraction',
                   'sb'    :'/Fluid Fields/Auxiliary/Entropy Per Baryon',
                   'xHeavy':'/Fluid Fields/Auxiliary/Heavy Mass Fraction',
                   'mu_n'  :'/Fluid Fields/Auxiliary/Neutron Chemical Potential',
                   'xn'    :'/Fluid Fields/Auxiliary/Neutron Mass Fraction',
                   'xp'    :'/Fluid Fields/Auxiliary/Proton Mass Fraction',
                   'press' :'/Fluid Fields/Auxiliary/Pressure',
                   'mu_p'  :'/Fluid Fields/Auxiliary/Proton Chemical Potential',
                   'gamma' :'/Fluid Fields/Auxiliary/Ratio of Specific Heats (Gamma)',
                   'c_s'   :'/Fluid Fields/Auxiliary/Sound Speed',
                   'u'     :'/Fluid Fields/Auxiliary/Specific Internal Energy',
                   'T'     :'/Fluid Fields/Auxiliary/Temperature',
                   'maxE'  :'/Fluid Fields/Diagnostic/Max E',
                   'minE'  :'/Fluid Fields/Diagnostic/Min E',
                   'shock1':'/Fluid Fields/Diagnostic/Shock (X1)',
                   'shock2':'/Fluid Fields/Diagnostic/Shock (X2)',
                   'shock3':'/Fluid Fields/Diagnostic/Shock (X3)',
                   'TCI'   :'/Fluid Fields/Diagnostic/TCI',
                   'theta1':'/Fluid Fields/Diagnostic/Theta 1',
                   'theta2':'/Fluid Fields/Diagnostic/Theta 2',
                   'theta3':'/Fluid Fields/Diagnostic/Theta 3'
                  }

    rad_dsets   = {
                   'N_e'   :'/Radiation Fields/Species_01/Conserved/Eulerian Number Density',
                   'N_eb'  :'/Radiation Fields/Species_02/Conserved/Eulerian Number Density',
                   'N_m'   :'/Radiation Fields/Species_03/Conserved/Eulerian Number Density',
                   'N_mb'  :'/Radiation Fields/Species_04/Conserved/Eulerian Number Density',
                   'N_t'   :'/Radiation Fields/Species_05/Conserved/Eulerian Number Density',
                   'N_tb'  :'/Radiation Fields/Species_06/Conserved/Eulerian Number Density',
                   'G_e1'  :'/Radiation Fields/Species_01/Conserved/Eulerian Number Flux Density (1)',
                   'G_e2'  :'/Radiation Fields/Species_01/Conserved/Eulerian Number Flux Density (2)',
                   'G_e3'  :'/Radiation Fields/Species_01/Conserved/Eulerian Number Flux Density (3)',
                   'G_eb1' :'/Radiation Fields/Species_02/Conserved/Eulerian Number Flux Density (1)',
                   'G_eb2' :'/Radiation Fields/Species_02/Conserved/Eulerian Number Flux Density (2)',
                   'G_eb3' :'/Radiation Fields/Species_02/Conserved/Eulerian Number Flux Density (3)',
                   'G_m1'  :'/Radiation Fields/Species_03/Conserved/Eulerian Number Flux Density (1)',
                   'G_m2'  :'/Radiation Fields/Species_03/Conserved/Eulerian Number Flux Density (2)',
                   'G_m3'  :'/Radiation Fields/Species_03/Conserved/Eulerian Number Flux Density (3)',
                   'G_mb1' :'/Radiation Fields/Species_04/Conserved/Eulerian Number Flux Density (1)',
                   'G_mb2' :'/Radiation Fields/Species_04/Conserved/Eulerian Number Flux Density (2)',
                   'G_mb3' :'/Radiation Fields/Species_04/Conserved/Eulerian Number Flux Density (3)',
                   'G_t1'  :'/Radiation Fields/Species_05/Conserved/Eulerian Number Flux Density (1)',
                   'G_t2'  :'/Radiation Fields/Species_05/Conserved/Eulerian Number Flux Density (2)',
                   'G_t3'  :'/Radiation Fields/Species_05/Conserved/Eulerian Number Flux Density (3)',
                   'G_tb1' :'/Radiation Fields/Species_06/Conserved/Eulerian Number Flux Density (1)',
                   'G_tb2' :'/Radiation Fields/Species_06/Conserved/Eulerian Number Flux Density (2)',
                   'G_tb3' :'/Radiation Fields/Species_06/Conserved/Eulerian Number Flux Density (3)',
                   'D_e'   :'/Radiation Fields/Species_01/Primitive/Lagrangian Number Density',
                   'D_eb'  :'/Radiation Fields/Species_02/Primitive/Lagrangian Number Density',
                   'D_m'   :'/Radiation Fields/Species_03/Primitive/Lagrangian Number Density',
                   'D_mb'  :'/Radiation Fields/Species_04/Primitive/Lagrangian Number Density',
                   'D_t'   :'/Radiation Fields/Species_05/Primitive/Lagrangian Number Density',
                   'D_tb'  :'/Radiation Fields/Species_06/Primitive/Lagrangian Number Density',
                   'I_e1'  :'/Radiation Fields/Species_01/Primitive/Lagrangian Number Flux Density (1)',
                   'I_e2'  :'/Radiation Fields/Species_01/Primitive/Lagrangian Number Flux Density (2)',
                   'I_e3'  :'/Radiation Fields/Species_01/Primitive/Lagrangian Number Flux Density (3)',
                   'I_eb1' :'/Radiation Fields/Species_02/Primitive/Lagrangian Number Flux Density (1)',
                   'I_eb2' :'/Radiation Fields/Species_02/Primitive/Lagrangian Number Flux Density (2)',
                   'I_eb3' :'/Radiation Fields/Species_02/Primitive/Lagrangian Number Flux Density (3)',
                   'I_m1'  :'/Radiation Fields/Species_03/Primitive/Lagrangian Number Flux Density (1)',
                   'I_m2'  :'/Radiation Fields/Species_03/Primitive/Lagrangian Number Flux Density (2)',
                   'I_m3'  :'/Radiation Fields/Species_03/Primitive/Lagrangian Number Flux Density (3)',
                   'I_mb1' :'/Radiation Fields/Species_04/Primitive/Lagrangian Number Flux Density (1)',
                   'I_mb2' :'/Radiation Fields/Species_04/Primitive/Lagrangian Number Flux Density (2)',
                   'I_mb3' :'/Radiation Fields/Species_04/Primitive/Lagrangian Number Flux Density (3)',
                   'I_t1'  :'/Radiation Fields/Species_05/Primitive/Lagrangian Number Flux Density (1)',
                   'I_t2'  :'/Radiation Fields/Species_05/Primitive/Lagrangian Number Flux Density (2)',
                   'I_t3'  :'/Radiation Fields/Species_05/Primitive/Lagrangian Number Flux Density (3)',
                   'I_tb1' :'/Radiation Fields/Species_06/Primitive/Lagrangian Number Flux Density (1)',
                   'I_tb2' :'/Radiation Fields/Species_06/Primitive/Lagrangian Number Flux Density (2)',
                   'I_tb3' :'/Radiation Fields/Species_06/Primitive/Lagrangian Number Flux Density (3)',
                   'Psi_e' :'/Radiation Fields/Species_01/Auxiliary/Lagrangian Eddington Factor',
                   'Psi_eb':'/Radiation Fields/Species_02/Auxiliary/Lagrangian Eddington Factor',
                   'Psi_m' :'/Radiation Fields/Species_03/Auxiliary/Lagrangian Eddington Factor',
                   'Psi_mb':'/Radiation Fields/Species_04/Auxiliary/Lagrangian Eddington Factor',
                   'Psi_t' :'/Radiation Fields/Species_05/Auxiliary/Lagrangian Eddington Factor',
                   'Psi_tb':'/Radiation Fields/Species_06/Auxiliary/Lagrangian Eddington Factor',
                   'h_e'   :'/Radiation Fields/Species_01/Auxiliary/Lagrangian Flux Factor',
                   'h_eb'  :'/Radiation Fields/Species_02/Auxiliary/Lagrangian Flux Factor',
                   'h_m'   :'/Radiation Fields/Species_03/Auxiliary/Lagrangian Flux Factor',
                   'h_mb'  :'/Radiation Fields/Species_04/Auxiliary/Lagrangian Flux Factor',
                   'h_t'   :'/Radiation Fields/Species_05/Auxiliary/Lagrangian Flux Factor',
                   'h_tb'  :'/Radiation Fields/Species_06/Auxiliary/Lagrangian Flux Factor',
                   'In_it' :'/Radiation Fields/Diagnostic/Inner Iterations',
                   'Out_it':'/Radiation Fields/Diagnostic/Outer Iterations',
                   'PL_dE' :'/Radiation Fields/Diagnostic/Positivity Limiter Energy Change',
                   'PL_th1':'/Radiation Fields/Diagnostic/Positivity Limiter Theta 1',
                   'PL_th2':'/Radiation Fields/Diagnostic/Positivity Limiter Theta 2'
                   }

    geom_dsets  = {
                   'cf'   :'/Geometry Fields/Conformal Factor',
                   'alp'  :'/Geometry Fields/Lapse Function',
                   'beta1':'/Geometry Fields/Shift Vector (1)',
                   'beta2':'/Geometry Fields/Shift Vector (2)',
                   'beta3':'/Geometry Fields/Shift Vector (3)',
                   'g11'  :'/Geometry Fields/Spatial Metric Component (11)',
                   'g22'  :'/Geometry Fields/Spatial Metric Component (22)',
                   'g33'  :'/Geometry Fields/Spatial Metric Component (33)',
                   'K11'  :'/Geometry Fields/Extrinsic Curvature Comp. (11)',
                   'K12'  :'/Geometry Fields/Extrinsic Curvature Comp. (12)',
                   'K13'  :'/Geometry Fields/Extrinsic Curvature Comp. (13)',
                   'K22'  :'/Geometry Fields/Extrinsic Curvature Comp. (22)',
                   'K23'  :'/Geometry Fields/Extrinsic Curvature Comp. (23)',
                   'K33'  :'/Geometry Fields/Extrinsic Curvature Comp. (33)',
                   'sdetg':'/Geometry Fields/Sqrt Spatial Metric Determinant',
                   'h1'   :'/Geometry Fields/Spatial Scale Factor (1)',
                   'h2'   :'/Geometry Fields/Spatial Scale Factor (2)',
                   'h3'   :'/Geometry Fields/Spatial Scale Factor (3)',
                   'Phi'  :'/Geometry Fields/Newtonian Potential'
                  }

    filelist = []
    keys     = []

    geom_filelist  = []
    geom_keys      = []
    fluid_filelist = []
    fluid_keys     = []
    rad_filelist   = []
    rad_keys       = []

    have_geom_data  = False
    have_fluid_data = False
    have_rad_data   = False

    n_files  = 0

    #time vector containing times from hdf5 output data files
    time     = []

    #energy total nodes and elements
    nodesE  = []
    cellsE   = []
    #total number of energy elements
    nE       = 0
    #number of energy nodes per element
    nNodesE  = 0

    #spatial total nodes and elements for each dimension
    nodesX1 = []
    cellsX1  = []
    nodesX2 = []
    cellsX2  = []
    nodesX3 = []
    cellsX3  = []
    #number of spatial elements in each dimension
    nX1 = 0; nX2 = 0; nX3 = 0
    #number of spatial nodes per element in each dimension
    nNodesX1 = 0; nNodesX2 = 0; nNodesX3 = 0

    #node number tables for
    #Radiation (4,nDOFZ), Z={E,X1,X2,X3}
    #and
    #Fluid     (3,nDOFX), X={X1,X2,X3}
    #fields
    NodeNumberTableX = []
    NodeNumberTableZ = []

    #Legendre-Gauss nodes and weights
    xLG_E  = []
    wLG_E  = []
    xLG_X1 = []
    wLG_X1 = []
    xLG_X2 = []
    wLG_X2 = []
    xLG_X3 = []
    wLG_X3 = []


    #weights to integrate elements
    WeightsX_q = []
    WeightsZ_q = []

    #number of degrees of freedom for energy(E),
    #spatial(X) and spatial x energy(Z)
    nDOFE = 0
    nDOFX = 0
    nDOFZ = 0

    def __init__(self, path=None, **kwargs):
        if path is None:
            raise ValueError("""Need to supply path to hdf5 datasets""")

        self.path     = path
        #store ordered lists of all geom, fluid, and radiation hdf5
        #files if present
        self.get_filelists(self.path)

        self.have_geom_data  = len(self.geom_filelist)
        self.have_fluid_data = len(self.fluid_filelist)
        self.have_rad_data   = len(self.rad_filelist)

        if not self.have_geom_data and not self.have_fluid_data and not self.have_rad_data:
            raise ValueError("""There are no hdf5 datasets
                                in the selected directory!!""")

        #open iteration 0 files and store all available keys
        if self.have_geom_data:
            f0 = h5.File(self.geom_filelist[0],'r')
            f0.visit(self.geom_keys.append)
            f0.close()
        if self.have_fluid_data:
            f0 = h5.File(self.fluid_filelist[0],'r')
            f0.visit(self.fluid_keys.append)
            f0.close()
        if self.have_rad_data:
            f0 = h5.File(self.rad_filelist[0],'r')
            f0.visit(self.rad_keys.append)
            f0.close()

        for key in self.time_dsets.keys():
            self.get_dset(key,0)
        for key in self.energy_dsets.keys():
            self.get_dset(key,0)
        for key in self.grid_dsets.keys():
            self.get_dset(key,0)
        for key in self.geom_dsets.keys():
            self.get_dset(key,0)
        for key in self.fluid_dsets.keys():
            self.get_dset(key,0)
        for key in self.rad_dsets.keys():
            self.get_dset(key,0)

        #if we don't actually want all data files, but
        #a maximum number of files
        if 'max_files' in kwargs:
            self.n_files = kwargs['max_files']
        else:
            if   self.have_geom_data:
                self.n_files = len(self.geom_filelist)
            elif self.have_fluid_data:
                self.n_files = len(self.fluid_filelist)
            elif self.have_rad_data:
                self.n_files = len(self.rad_filelist)

        if self.have_geom_data and self.have_fluid_data \
        and len(self.geom_filelist) != len(self.fluid_filelist):
            raise ValueError("""geom and fluid datasets do not
                                have the same number of files""")
        if self.have_geom_data and self.have_rad_data \
        and len(self.geom_filelist) != len(self.rad_filelist):
            raise ValueError("""geom and radiation datasets do not
                                have the same number of files""")
        if self.have_fluid_data and self.have_rad_data \
        and len(self.fluid_filelist) != len(self.rad_filelist):
            raise ValueError("""fluid and radiation datasets do not
                                have the same number of files""")
        if self.have_geom_data and self.have_fluid_data and self.have_rad_data \
        and not len(self.geom_filelist) == len(self.fluid_filelist) == len(self.rad_filelist):
            raise ValueError("""geom, fluid, and radiation  datasets do not
                                have the same number of files""")

        self.calc_nNodesX_and_nDOFX()
        self.calc_NodeNumberTableX()
        self.calc_WeightsX_q()

        if self.have_rad_data:
            self.calc_nNodes_and_nDOFZ()
            self.calc_NodeNumberTableZ()
            self.calc_WeightsZ_q()

    #get a list of files from path+filename wildcards
    def get_filelists(self, path):
        self.geom_filelist  = sorted(glob.glob(path+'*GeometryFields_*.h5'))
        self.fluid_filelist = sorted(glob.glob(path+'*FluidFields_*.h5'))
        self.rad_filelist   = sorted(glob.glob(path+'*RadiationFields_*.h5'))

    #load dataset into numpy array at given filenumer/iteration
    #provide the name of the wanted dataset which has to be
    #in one of the following class dictionaries:
    #time_dsets, grid_dsets, energy_dsets,
    #geom_dsets, fluid_dsets, rad_dsets
    #we then go through a series of if
    #statements to check where to get the
    #dataset from, as some data (such as time and
    #the spatial grid) are present in all field
    #type output files
    def get_dset(self, dataset, iteration):

        if   dataset in self.geom_dsets.keys():
            f = h5.File(self.geom_filelist[iteration], 'r')
            dset = np.asarray(f[self.geom_dsets[dataset]])
            f.close
        elif dataset in self.fluid_dsets.keys():
            f = h5.File(self.fluid_filelist[iteration], 'r')
            dset = np.asarray(f[self.fluid_dsets[dataset]])
            f.close
        elif dataset in self.rad_dsets.keys():
            f = h5.File(self.rad_filelist[iteration], 'r')
            dset = np.asarray(f[self.rad_dsets[dataset]])
            f.close
        #energy grid exists only in Radiation files
        elif dataset in self.energy_dsets.keys():
            f = h5.File(self.rad_filelist[iteration], 'r')
            dset = np.asarray(f[self.energy_dsets[dataset]])
            f.close
        elif dataset in self.time_dsets.keys():
            if   self.have_geom_data:
                f = h5.File(self.geom_filelist[iteration], 'r')
                dset = np.asarray(f[self.time_dsets[dataset]])
                f.close
            elif self.have_fluid_data:
                f = h5.File(self.fluid_filelist[iteration], 'r')
                dset = np.asarray(f[self.time_dsets[dataset]])
                f.close
            elif self.have_rad_data:
                f = h5.File(self.rad_filelist[iteration], 'r')
                dset = np.asarray(f[self.time_dsets[dataset]])
                f.close
        elif dataset in self.grid_dsets.keys():
            if   self.have_geom_data:
                f = h5.File(self.geom_filelist[iteration], 'r')
                dset = np.asarray(f[self.grid_dsets[dataset]])
                f.close
            elif self.have_fluid_data:
                f = h5.File(self.fluid_filelist[iteration], 'r')
                dset = np.asarray(f[self.grid_dsets[dataset]])
                f.close
            elif self.have_rad_data:
                f = h5.File(self.rad_filelist[iteration], 'r')
                dset = np.asarray(f[self.grid_dsets[dataset]])
                f.close
        #dataset wasn't found...
        else:
            print('dataset ',dataset, ' was not found')
            raise ValueError("""Error in get_dset!""")
        return dset

    #show all available datasets in hdf5 files
    def show_available_datasets(self, fields):
        if fields == 'Geometry':
            print(self.geom_keys)
        elif fields == 'Fluid':
            print(self.fluid_keys)
        elif fields == 'Radiation':
            print(self.rad_keys)
        else:
            raise ValueError("""Unsupported field type
                             in show_available_datasets!""")

    #store snapshots in time as +1 dimensional numpy array
    def get_timeseries(self, dataset=None):
        timeseries = []
        for it in range(self.n_files):
            data = self.get_dset(dataset, it)
            timeseries.append(data)
        return np.asarray(timeseries)

    #get timevector from datafiles
    def get_time_vector(self):
        time = []
        for it in range(self.n_files):
            data = self.get_dset('t', it)
            time.append(data)
        self.time = np.asarray(time)
        self.time = self.time.ravel()
        return self.time

    #get energy nodes
    def get_energy_nodes(self):
       if self.have_rad_data:
           self.nodesE = self.get_dset('E', 0)
       else:
           raise ValueError("""No Radiation fields data
                               files found to get energy nodes""")

    #get energy grid
    def get_energy_grid(self):
       if self.have_rad_data:
           self.cellsE = self.get_dset('E_C',0)
       else:
           raise ValueError("""No Radiation fields data
                               files found to get energy grid""")

    #get spatial grid nodes
    def get_spatial_nodes(self):
        self.nodesX1 = self.get_dset('X1', 0)
        self.nodesX2 = self.get_dset('X2', 0)
        self.nodesX3 = self.get_dset('X3', 0)

    #get spatial grid
    def get_spatial_grid(self):
        self.cellsX1 = self.get_dset('X1_C',0)
        self.cellsX2 = self.get_dset('X2_C',0)
        self.cellsX3 = self.get_dset('X3_C',0)

    def calc_nNodesX_and_nDOFX(self):

        self.get_spatial_nodes()
        self.get_spatial_grid()

        self.nX1 = self.cellsX1.size
        self.nX2 = self.cellsX2.size
        self.nX3 = self.cellsX3.size

        self.nNodesX1 = self.nodesX1.size//self.nX1
        self.nNodesX2 = self.nodesX2.size//self.nX2
        self.nNodesX3 = self.nodesX3.size//self.nX3

        self.nDOFX = self.nNodesX1 \
                   * self.nNodesX2 \
                   * self.nNodesX3

    def calc_nNodes_and_nDOFZ(self):

        self.get_energy_nodes()
        self.get_energy_grid()

        self.nE      = self.cellsE.size
        self.nNodesE = self.nodesE.size//self.nE

        self.nDOFE = self.nNodesE
        self.nDOFX = self.nNodesX1 \
                   * self.nNodesX2 \
                   * self.nNodesX3
        self.nDOFZ = self.nNodesX1 \
                   * self.nNodesX2 \
                   * self.nNodesX3 \
                   * self.nNodesE

    def calc_NodeNumberTableX(self):
        iNode = 0
        NodeNumberTable = np.full((3,self.nDOFX),np.nan)
        for iNodeX3 in range(1,self.nNodesX3+1):
            for iNodeX2 in range(1,self.nNodesX2+1):
                for iNodeX1 in range(1,self.nNodesX1+1):
                    iNode += 1
                    NodeNumberTable[:,iNode-1] = [iNodeX1-1, iNodeX2-1, iNodeX3-1]

        if np.sum(NodeNumberTable) == np.nan:
            raise ValueError("""did not calculate NodeNumberTableX array correctly""")
        self.NodeNumberTableX = NodeNumberTable.astype(int)

    def calc_NodeNumberTableZ(self):
        iNode = 0
        NodeNumberTable = np.full((4,self.nDOFZ),np.nan)
        for iNodeX3 in range(1,self.nNodesX3+1):
            for iNodeX2 in range(1,self.nNodesX2+1):
                for iNodeX1 in range(1,self.nNodesX1+1):
                    for iNodeE in range(1,self.nNodesE+1):
                        iNode += 1
                        NodeNumberTable[:,iNode-1] = [iNodeE-1, iNodeX1-1, iNodeX2-1, iNodeX3-1]

        if np.sum(NodeNumberTable) == np.nan:
            raise ValueError("""did not calculate NodeNumberTableZ array correctly""")
        self.NodeNumberTableZ = NodeNumberTable.astype(int)

    def calc_WeightsX_q(self):
        self.xLG_X1, self.wLG_X1 = self.get_quadrature(self.nNodesX1)
        self.xLG_X2, self.wLG_X2 = self.get_quadrature(self.nNodesX2)
        self.xLG_X3, self.wLG_X3 = self.get_quadrature(self.nNodesX3)

        weights_q = np.full((self.nDOFX), np.nan)

        for iNode in range(1,self.nDOFX+1):
            iNodeX1 = self.NodeNumberTableX[0,iNode-1]
            iNodeX2 = self.NodeNumberTableX[1,iNode-1]
            iNodeX3 = self.NodeNumberTableX[2,iNode-1]

            weights_q[iNode-1] = self.wLG_X1[iNodeX1-1] \
                               * self.wLG_X2[iNodeX2-1] \
                               * self.wLG_X3[iNodeX3-1]

        if np.sum(weights_q) == np.nan:
            raise ValueError("""did not calculate WeightsX_q array correctly""")

        self.WeightsX_q = weights_q

    def calc_WeightsZ_q(self):
        self.xLG_X1, self.wLG_X1 = self.get_quadrature(self.nNodesX1)
        self.xLG_X2, self.wLG_X2 = self.get_quadrature(self.nNodesX2)
        self.xLG_X3, self.wLG_X3 = self.get_quadrature(self.nNodesX3)

        self.xLG_E, self.wLG_E = self.get_quadrature(self.nNodesE)

        weights_q = np.full((self.nDOFZ), np.nan)

        for iNode in range(1,self.nDOFZ+1):
            iNodeE  = self.NodeNumberTableZ[0,iNode-1]
            iNodeX1 = self.NodeNumberTableZ[1,iNode-1]
            iNodeX2 = self.NodeNumberTableZ[2,iNode-1]
            iNodeX3 = self.NodeNumberTableZ[3,iNode-1]

            weights_q[iNode-1] = self.wLG_E[iNodeE-1]   \
                               * self.wLG_X1[iNodeX1-1] \
                               * self.wLG_X2[iNodeX2-1] \
                               * self.wLG_X3[iNodeX3-1]
        if np.sum(weights_q) == np.nan:
            raise ValueError("""did not calculate WeightsZ_q array correctly""")

        self.WeightsZ_q = weights_q

    def get_quadrature(self, nNodes):
        xNodes, wNodes = np.polynomial.legendre.leggauss(nNodes)
        xNodes /= 2.0; wNodes /=2.0
        return xNodes, wNodes

    def reshape_dataE(self, dset):
        data = np.full((self.nDOFE, self.nE), np.nan)

        for iE in range(1,self.nE+1):
            for iNodeE in range(1,self.nDOFE+1):
                iN_E   = iNodeE  + (iE -1)*self.nNodesE
                data[iNodeE-1,iE-1] = dset[iN_E-1]
        if np.sum(data) == np.nan:
            raise ValueError("""did not reshape_dataE correctly""")
        return data

    #the Fluid and Geometry field data in the hdf5 files is stored as
    #(nNodesX3*nX3, nNodesX2*nX2, nNodesX1*nX1)
    #in order to use existing thornado routines to integrate
    #over nodes in elements, we need to reshape this data into
    #the way it is stored in thornado:
    #(nDOFX, nX1, nX2, nX3)
    def reshape_dataX(self, dset):
        data = np.full((self.nDOFX, self.nX1, self.nX2, self.nX3), np.nan)

        for iX3 in range(1,self.nX3+1):
            for iX2 in range(1,self.nX2+1):
                for iX1 in range(1,self.nX1+1):
                    for iNodeX3 in range(1,self.nNodesX3+1):
                        for iNodeX2 in range(1,self.nNodesX2+1):
                            for iNodeX1 in range(1,self.nNodesX1+1):
                                iNodeX = iNodeX1 + (iNodeX2-1)*self.nNodesX1 \
                                       + (iNodeX3-1)*self.nNodesX1*self.nNodesX2
                                iN_X1  = iNodeX1 + (iX1-1)*self.nNodesX1
                                iN_X2  = iNodeX2 + (iX2-1)*self.nNodesX2
                                iN_X3  = iNodeX3 + (iX3-1)*self.nNodesX3
                                data[iNodeX-1,iX1-1,iX2-1,iX3-1] = dset[iN_X3-1,iN_X2-1,iN_X1-1]

        if np.sum(data) == np.nan:
            raise ValueError("""did not reshape_dataX correctly""")
        return data

    #the Radiation field data in the hdf5 files is stored as
    #(nNodesX3*nX3, nNodesX2*nX2, nNodesX1*nX1, nNodesE*nE)
    #in order to use existing thornado routines to integrate
    #over nodes in elements, we need to reshape this data into
    #the way it is stored in thornado:
    #(nDOFZ,nE, nX1, nX2, nX3)
    def reshape_dataZ(self, dset):
        data = np.full((self.nDOFZ,self.nE,self.nX1,self.nX2,self.nX3),np.nan)

        for iX3 in range(1,self.nX3+1):
            for iX2 in range(1,self.nX2+1):
                for iX1 in range(1,self.nX1+1):
                    for iE in range(1,self.nE+1):
                        for iNodeX3 in range(1,self.nNodesX3+1):
                            for iNodeX2 in range(1,self.nNodesX2+1):
                                for iNodeX1 in range(1,self.nNodesX1+1):
                                    for iNodeE in range(1,self.nNodesE+1):
                                        iNodeZ = iNodeE + (iNodeX1-1)*self.nNodesE \
                                               + (iNodeX2-1)*self.nNodesE*self.nNodesX1 \
                                               + (iNodeX3-1)*self.nNodesE*self.nNodesX1*self.nNodesX2
                                        iN_E   = iNodeE  + (iE -1)*self.nNodesE
                                        iN_X1  = iNodeX1 + (iX1-1)*self.nNodesX1
                                        iN_X2  = iNodeX2 + (iX2-1)*self.nNodesX2
                                        iN_X3  = iNodeX3 + (iX3-1)*self.nNodesX3
                                        data[iNodeZ-1,iE-1,iX1-1,iX2-1,iX3-1] = dset[iN_X3-1,iN_X2-1,iN_X1-1,iN_E-1]

        if np.sum(data) == np.nan:
            raise ValueError("""did not reshape_dataZ correctly""")
        return data

    #get spatial element averages from thornado spatial nodal points
    def get_element_averageX(self, dataset, iteration, use_metric=None):

        h5_data = self.get_dset(dataset, iteration)

        data = self.reshape_dataX(h5_data)

        element_average = np.full((self.nX1, self.nX2, self.nX3),np.nan)

        if use_metric:

            if not self.have_geom_data:
                raise ValueError("""Requested sqrt of metric determiant
                                   to calculate element_averageZ, but
                                   no geometry fields hdf5 files were found!""")

            h5_sdetg = self.get_dset('sdetg', iteration)
            sdetg = self.reshape_dataX(h5_sdetg)

            for iX3 in range(1,self.nX3+1):
                for iX2 in range(1,self.nX2+1):
                    for iX1 in range(1,self.nX1+1):
                        sum1 = 0.0
                        sum2 = 0.0
                        for iNodeX in range(1,self.nDOFX+1):
                            sum1 += self.WeightsX_q[iNodeX-1] \
                                  * data[iNodeX-1,iX1-1,iX2-1,iX3-1] \
                                  * sdetg[iNodeX-1,iX1-1,iX2-1,iX3-1]
                            sum2 += self.WeightsX_q[iNodeX-1] \
                                  * sdetg[iNodeX-1,iX1-1,iX2-1,iX3-1]
                        element_average[iX1-1,iX2-1,iX3-1] = sum1/sum2

            if np.sum(element_average) == np.nan:
                raise ValueError("""did not calculate
                                    element_averageX array correctly""")

            return element_average

        else:

            for iX3 in range(1,self.nX3+1):
                for iX2 in range(1,self.nX2+1):
                    for iX1 in range(1,self.nX1+1):
                        sum = 0.0
                        for iNodeX in range(1,self.nDOFX+1):
                            sum += self.WeightsX_q[iNodeX-1]*data[iNodeX-1,iX1-1,iX2-1,iX3-1]
                        element_average[iX1-1,iX2-1,iX3-1] = sum

            if np.sum(element_average) == np.nan:
                raise ValueError("""did not calculate
                                    element_averageX array correctly""")

            return element_average

    #get energy cell averages from thornado energy and spatial nodal points
    def get_element_averageZ(self, dataset, iteration, use_metric=None):

        h5_data = self.get_dset(dataset, iteration)

        data = self.reshape_dataZ(h5_data)

        element_average = np.full((self.nE, self.nX1, self.nX2, self.nX3),np.nan)

        if use_metric:

            if not self.have_geom_data:
                raise ValueError("""Requested sqrt of metric determiant
                                   to calculate element_averageZ, but
                                   no geometry fields hdf5 files were found!""")

            h5_sdetg = self.get_dset('sdetg', iteration)
            sdetg    = self.reshape_dataX(h5_sdetg)
            E        = self.reshape_dataE(self.nodesE)
            Tau_q    = np.full((self.nDOFZ),np.nan)

            for iX3 in range(1,self.nX3+1):
                for iX2 in range(1,self.nX2+1):
                    for iX1 in range(1,self.nX1+1):
                        for iE in range(1,self.nE+1):
                            for iNodeX in range(1,self.nDOFX+1):
                                for iNodeE in range(1,self.nDOFE+1):
                                    iNode = iNodeE+(iNodeX-1)*self.nDOFE
                                    Tau_q[iNode-1] = sdetg[iNodeX-1,iX1-1,iX2-1,iX3-1] \
                                                   * E[iNodeE-1,iE-1]**2
                            sum1 = 0.0
                            sum2 = 0.0
                            for iNodeZ in range(1,self.nDOFZ+1):
                                sum1 += self.WeightsZ_q[iNodeZ-1]*data[iNodeZ-1,iE-1,iX1-1,iX2-1,iX3-1] \
                                      * Tau_q[iNodeZ-1]
                                sum2 += self.WeightsZ_q[iNodeZ-1]*Tau_q[iNodeZ-1]
                            element_average[iE-1,iX1-1,iX2-1,iX3-1] = sum1/sum2

            if np.sum(element_average) == np.nan:
                raise ValueError("""did not calculate
                                    element_averageZ array correctly""")

            return element_average

        else:

            for iX3 in range(1,self.nX3+1):
                for iX2 in range(1,self.nX2+1):
                    for iX1 in range(1,self.nX1+1):
                        for iE in range(1,self.nE+1):
                            sum = 0.0
                            for iNodeZ in range(1,self.nDOFZ+1):
                                sum += self.WeightsZ_q[iNodeZ-1]*data[iNodeZ-1,iE-1,iX1-1,iX2-1,iX3-1]
                            element_average[iE-1,iX1-1,iX2-1,iX3-1] = sum

            if np.sum(element_average) == np.nan:
                raise ValueError("""did not calculate
                                    element_averageZ array correctly""")

            return element_average

#example use

#! from thornady_hdf5 import thornado_hdf5 as h5

#data = h5(path='path_to_files')
#nodesE = data.nodesE
#cellsE = data.cellsE

#t = data.get_time_vector()

#its = t.size
#D_e                 = data.get_dset('D_e',its-1)
#cell_averages_D_e   = data.get_element_averageZ('D_e',its-1)
#cell_avg_metric_D_e = data.get_element_averageZ('D_e',its-1,use_metric=True)
