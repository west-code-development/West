##
## Copyright (C) 2015 M. Govoni 
## This file is distributed under the terms of the
## GNU General Public License. See the file `License'
## in the root directory of the present distribution,
## or http://www.gnu.org/copyleft/gpl.txt .
##
## This file is part of WEST
##
## Contributors to this file: 
## Marco Govoni
##

from __future__ import print_function

##
## Header 
##

def header() : 
   import datetime
   print(" _    _ _____ _____ _____            ")       
   print("| |  | |  ___/  ___|_   _|           ")
   print("| |  | | |__ \ `--.  | |_ __  _   _  ")
   print("| |/\| |  __| `--. \ | | '_ \| | | | ")
   print("\  /\  / |___/\__/ / | | |_) | |_| | ")
   print(" \/  \/\____/\____/  \_/ .__/ \__, | ")
   print("                       | |     __/ | ")
   print("                       |_|    |___/  ")
   print("*** <WEST> ***")
   print("WEST version     : 3.0.0")
   print("Today            : ", datetime.datetime.today())
   print("*** </WEST> ***")
   print(" ")

##
## Versions 
##

def versions() :
   import numpy as np
   import scipy as sp
   import matplotlib as mpl
   print("*** <VERSIONS> ***")
   print("Numpy version      :", np.__version__)
   print("Scipy version      :", sp.__version__)
   print("MatPlotLib version :", mpl.__version__)
   print("*** </VERSION> ***")
   print(" ")

##
## Footer
##

def footer() : 
   import datetime
   print("*** <END> ***")
   print("Job done : ", datetime.datetime.today())
   print("*** </END> ***")

##
## Constants
##

atomic_symbol = ( "blank", "H", "He", "Li", "Be", "B", "C", "N", "O", "F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db" ,"Sg" ,"Bh" ,"Hs" ,"Mt" ,"Ds","Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og" )

##
## Useful functions 
##

def gaussian(x, mu, sig):
    import numpy as np
    return 1./(np.sqrt(2.*np.pi)*sig)*np.exp(-np.power((x - mu)/sig, 2.)/2)

##
## CubeFile
##

class CubeFile:
    #
    def __init__(self, name):
        self.fname = name
        self.nat = 0 
        self.shift = []
        self.n = []
        self.acell = []
        self.atsymbol = []
        self.atcharge = []
        self.atpos = []
        self.data = []
        self.nspec = 0
        self.spec_list = [] 
        self.read() 
    #
    def read(self):
        import numpy as np
        print('Reading Cube File : ', self.fname) 
        # read
        file_handle = open(self.fname, 'r')
        lines_list = file_handle.readlines()
        nline = 0 
        # skip first 2 lines 
        nline += 2
        # nat, shift (1 line)
        my_line = lines_list[nline].split()
        self.shift = np.zeros( ( 3 ) )
        self.nat = int( my_line[0] )
        self.shift = [float(val) for val in my_line[1:]]
        nline += 1
        # n acell (3 lines)
        self.n = np.zeros( ( 3 ) )
        self.acell = np.zeros( ( 3, 3 ) )
        for i in range(3) : 
            my_line = lines_list[nline].split()
            self.n[i] = int( my_line[0] )
            self.acell[i,0:] = [float(val)*float(self.n[i]) for val in my_line[1:]]
            nline += 1 
        self.n=self.n.astype(int)
        # atomic allocations (nat lines) 
        self.atsymbol = np.zeros( ( self.nat ) )
        self.atcharge = np.zeros( ( self.nat ) )
        self.atpos = np.zeros( ( self.nat, 3 ) )
        for ia in range(self.nat) : 
            my_line = lines_list[nline].split()
            self.atsymbol[ia] = int( my_line[0] )
            self.atcharge[ia] = float( my_line[1] )
            self.atpos[ia,0:] = [float(val) for val in my_line[2:]]
            nline += 1 
        self.atsymbol=self.atsymbol.astype(int)
        # data allocations (xxx lines) 
        self.data = [[float(val) for val in line.split()] for line in lines_list[nline:]]
        self.data = np.reshape(self.data, (self.n[0],self.n[1],self.n[2]))
        #
        for iat in range( self.nat ) : 
            if (self.atsymbol[iat] not in self.spec_list) :
                self.spec_list.append(self.atsymbol[iat])
        self.nspec = len( self.spec_list )
    #
    def show_content(self):
        import numpy as np
        print('File : ', self.fname) 
        print(self.nat, self.shift[0], self.shift[1], self.shift[2])
        for i in range(3) : 
            print(self.n[i], self.acell[i,0:3])
        for ia in range(self.nat) : 
            print(atomic_symbol[int(self.atsymbol[ia])], self.atcharge[ia], self.atpos[ia,0:3])
        print(self.nspec, self.spec_list)
    #
    def interpolate(self,pts):
        import numpy as np
        data_periodic = np.zeros((self.n[0]+1,self.n[1]+1,self.n[2]+1))
        for i in range(self.n[0]+1):
            ii = i % self.n[0]
            for j in range(self.n[1]+1):
                jj = j % self.n[1]
                for k in range(self.n[2]+1):
                   kk = k % self.n[2]
                   data_periodic[i,j,k] = self.data[ii,jj,kk]
        lx = np.sqrt( np.dot(self.acell[0,:],self.acell[0,:]) ) 
        ly = np.sqrt( np.dot(self.acell[1,:],self.acell[1,:]) ) 
        lz = np.sqrt( np.dot(self.acell[2,:],self.acell[2,:]) ) 
        x = np.linspace(0, lx, self.n[0]+1, endpoint=True)
        y = np.linspace(0, ly, self.n[1]+1, endpoint=True)
        z = np.linspace(0, lz, self.n[2]+1, endpoint=True)
        #  
        from scipy.interpolate import RegularGridInterpolator
        my_interpolating_function = RegularGridInterpolator((x, y, z), data_periodic)
        #
        return my_interpolating_function(pts)
    #
    def get_slice(self,mesh2d,xmap,ymap,zmap):
        import numpy as np
        pts = np.zeros( ( (mesh2d[0])*(mesh2d[1]), 3 ) )
        t = 0
        for i in range(mesh2d[0]):
            for j in range(mesh2d[1]):
                pts[t,0] = xmap[i,j]
                pts[t,1] = ymap[i,j]
                pts[t,2] = zmap[i,j]  
                t += 1 
        cmap = self.interpolate(pts)
        cmap = np.reshape(cmap, (mesh2d[0],mesh2d[1]) )
        return cmap

##
## EnergyFile
##

class EnergyFile:
    #
    def __init__(self, name):
        self.fname = name
        self.ne = 0 
        self.e = []
        self.name = []
        self.has_names = False
        self.npt = 0
        self.plot = []
        self.axis = []
        self.read() 
    #
    def read(self):
        import numpy as np
        print('Reading Energy File : ', self.fname) 
        # read
        file_handle = open(self.fname, 'r')
        lines_list = file_handle.readlines()
        self.ne = 0
        for line in lines_list : 
            if not line.startswith('#'): 
               self.ne += 1
        self.e = np.zeros( ( self.ne, 3 ) )
        ie = 0 
        for line in lines_list :
            if not line.startswith('#'): 
               my_line = line.split() 
               self.e[ie,0] = float( my_line[0] ) 
               self.e[ie,1] = float( my_line[1] ) 
               self.e[ie,2] = float( my_line[2] )
               if( len(my_line) == 4 ) :
                  self.name.append( my_line[3] )
                  self.has_names = True
               ie+=1
    #
    def show_content(self):
        import numpy as np
        print('File : ', self.fname)
        print('Detected no. energies : ', self.ne) 
        if ( self.has_names ) :
           print('[Energy, Sigma, Weight], Name')
           for ie in range(self.ne) :
               print(self.e[ie,0:], self.name[ie])
        else :
           print('[Energy, Sigma, Weight]')
           for ie in range(self.ne) :
               print(self.e[ie,0:])
    #
    def crunch_datain(self):
        import numpy as np
        if ( self.has_names ) :
           file_handle = open(self.name[0], 'r')
           lines_list = file_handle.readlines()
           self.npt = 0
           for line in lines_list : 
               if not line.startswith('#'): 
                   self.npt += 1
           self.axis = np.zeros( ( self.npt ) )
           ipt = 0 
           for line in lines_list :
               if not line.startswith('#'): 
                  my_line = line.split()
                  self.axis[ipt] = float( my_line[0] )
                  ipt += 1
           self.plot = np.zeros( ( self.ne, self.npt ) )
           for ie in range(self.ne) :
               file_handle = open(self.name[ie], 'r')
               lines_list = file_handle.readlines()
               ipt = 0 
               for line in lines_list :
                   if not line.startswith('#'): 
                      my_line = line.split()
                      self.plot[ie,ipt] = float( my_line[1] )
                      ipt += 1 
        else :
           print('No data input')

##
## Plot DOS
##

def plot_dos(energy_file_name,output_file_name,erange) : 
    #
    ### USER DEFINITION ###
    #energy_file_name      = ['en1.dat', 'en2.dat']
    #output_file_name      = 'dos.png'
    #erange                = [40.,60.,0.01] eV
    ### END ##############
    #
    # N.B. 
    # the file energy_file_name should contain 3 values per line : en, sigma, weight 
    # erange[0] = emin 
    # erange[1] = emax 
    # erange[2] = deltae 
    #
    import numpy as np
    import scipy as sp
    import matplotlib as mpl
    #
    header()
    versions()
    #
    print("*** <INPUT> ***")
    print("Energy file name     : ", energy_file_name) 
    print("Output file name     : ", output_file_name)
    print("erange               : ", erange)
    print("*** </INPUT> ***")
    print(" ")
    # 
    print("*** <DOS> ***")
    #
    npte = int((erange[1]-erange[0])/erange[2]) + 1
    print("Using npte : ", npte)
    #
    e_axis = np.linspace(erange[0], erange[1], npte, endpoint=True)
    #
    dos_axis = {}
    emin = []
    emax = []
    ymax = []
    for ifile in energy_file_name : 
       dos_axis[ifile] = np.zeros ( (npte) )
       en = EnergyFile(ifile)
       en.show_content()
       #
       emin.append( np.min(en.e[:,0]) )
       emax.append( np.max(en.e[:,0]) )
       #
       for ie in range(en.ne) :
          for ix in range(npte) : 
             dos_axis[ifile][ix] += gaussian( e_axis[ix], en.e[ie,0], en.e[ie,1] ) * en.e[ie,2] 
       #
       ymax.append( np.max(dos_axis[ifile][:]) )
    #
    print("Requested (emin,emax) : ", erange[0], erange[1])
    print("Detected  (emin,emax) : ", np.min(emin), np.max(emax))
    #
    print("*** </DOS> ***")
    print(" ")
    #
    print("*** <PLOT> ***")
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    for ifile in energy_file_name : 
       dos_plot = ax.plot(e_axis,dos_axis[ifile],label=ifile)
       
    plt.xlim([erange[0],erange[1]])
    plt.ylim([0,np.max(ymax[:])])
    plt.xlabel('En (eV)')
    plt.ylabel('DOS')
    plt.savefig(output_file_name)
    plt.legend()
    print("output written in : ", output_file_name)
    print("waiting for user to close image preview...")
    plt.show()
    fig.clear()
    print("*** </PLOT> ***")
    print(" ")
    #
    footer()

##
## Plot 2D Slice
##

def plot_2dslice(cube_file_name,output_file_name,mesh2d,slice_direction,const,add_atoms) : 
    #  
    ### USER DEFINITION ###
    #cube_file_name   = 'band2.cub'
    #output_file_name = 'band2.jpeg' 
    #mesh2d           = [ 400, 400 ] 
    #slice_direction  = 2 # cn be 0, 1 or 2 
    #const            = 10    # in a.u.
    #add_atoms        = True
    ### END ##############
    #
    # N.B. 
    # The file cube_file_name should be a cube file. 
    # 
    import numpy as np
    import scipy as sp
    import matplotlib as mpl
    #
    header()
    versions()
    #
    print("*** <INPUT> ***")
    print("Cube file name   : ", cube_file_name)
    print("Output file name : ", output_file_name) 
    print("mesh2d           : ", mesh2d)
    print("Slice direction  : ", slice_direction)
    print("Cut (a.u.)       : ", const)
    print("add_atoms        : ", add_atoms) 
    print("*** </INPUT> ***")
    print(" ")
    #
    print("*** <SLICING> ***")
    print("using point (x,y) : ", mesh2d)
    #
    label = [ 'x', 'y', 'z' ]
    #
    if ( slice_direction == 0 ) : 
       dir_x = 1 
       dir_y = 2 
    if ( slice_direction == 1 ) : 
       dir_x = 0 
       dir_y = 2 
    if ( slice_direction == 2 ) : 
       dir_x = 0 
       dir_y = 1 
    #
    print("along direction   : ", label[slice_direction])
    #
    cube = CubeFile(cube_file_name)
    #
    directionx = cube.acell[dir_x,:] 
    directiony = cube.acell[dir_y,:]  
    directionc = cube.acell[slice_direction,:]
    lx = np.sqrt(np.dot(directionx,directionx))
    ly = np.sqrt(np.dot(directiony,directiony))
    lc = np.sqrt(np.dot(directionc,directionc))
    #
    xmap = np.zeros( ( (mesh2d[0]+1), (mesh2d[1]+1) ) )
    ymap = np.zeros( ( (mesh2d[0]+1), (mesh2d[1]+1) ) )
    zmap = np.zeros( ( (mesh2d[0]+1), (mesh2d[1]+1) ) )
    for i in range(mesh2d[0]+1):
        for j in range(mesh2d[1]+1):
            xmap[i,j] = float(i)/float(mesh2d[0])*lx
            ymap[i,j] = float(j)/float(mesh2d[1])*ly
            zmap[i,j] = const
    #
    cmap = cube.get_slice(mesh2d,xmap,ymap,zmap)
    #
    print("generated X,Y,C maps")
    print("*** </SLICING> ***")
    print(" ")
    #  
    print("*** <PLOT> ***")
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    heatmap = ax.pcolormesh(xmap,ymap,cmap)
    cbar = plt.colorbar(heatmap)
    plt.xlim([0,lx])
    plt.ylim([0,ly])
    plt.xlabel(label[dir_x] + ' (a.u.)')
    plt.ylabel(label[dir_y] + ' (a.u.)')
    plt.title(label[slice_direction] + " const = " + str(const) + " (a.u.)")
    plt.axis('image')
    if ( add_atoms ) : 
       for iat in range(cube.nat) :
           xau = cube.atpos[iat,dir_x] + cube.shift[dir_x] 
           if ( xau < 0. ) : 
              xau += lx 
           xau = xau % lx
           yau = cube.atpos[iat,dir_y] + cube.shift[dir_x] 
           if ( yau < 0. ) : 
              yau += ly 
           yau = yau % ly
           circ=plt.Circle((xau,yau), radius=max(lx,ly)/50, color='w', fill=True)
           plt.text(xau, yau, atomic_symbol[cube.atsymbol[iat]], fontsize = 11, color = 'r', horizontalalignment='center', verticalalignment='center')
           ax.add_patch(circ)
    plt.savefig(output_file_name)
    print("output written in : ", output_file_name)
    print("waiting for user to close image preview...")
    plt.show()
    fig.clear()
    print("*** </PLOT> ***")
    print(" ")
    #
    footer()

##
## Plot Abundacy of species 
##

def plot_species_abundancy(cube_file_name,output_file_name) : 
    #
    ### USER DEFINITION ###
    #cube_file_name        = 'band1.cube'
    #output_file_name      = 'spab.png'
    ### END ##############
    #
    import numpy as np
    import scipy as sp
    import matplotlib as mpl
    #
    header()
    versions()
    #
    print("*** <INPUT> ***")
    print("Cube file name       : ", cube_file_name) 
    print("Output file name     : ", output_file_name) 
    print("*** </INPUT> ***")
    print(" ")
    #
    print("*** <ABUNDANCIES> ***")
    #
    cube = CubeFile(cube_file_name)
    #
    lx = np.sqrt(np.dot(cube.acell[0,:],cube.acell[0,:]))
    ly = np.sqrt(np.dot(cube.acell[1,:],cube.acell[1,:]))
    lz = np.sqrt(np.dot(cube.acell[2,:],cube.acell[2,:]))
    # 
    x_axis = np.linspace(0, lx, 4*cube.n[0]+1, endpoint=True)
    y_axis = np.linspace(0, ly, 4*cube.n[1]+1, endpoint=True)
    z_axis = np.linspace(0, lz, 4*cube.n[2]+1, endpoint=True)
    #
    spec_x = np.zeros ( ( cube.nspec, 4*cube.n[0]+1 ) )
    spec_y = np.zeros ( ( cube.nspec, 4*cube.n[1]+1 ) )
    spec_z = np.zeros ( ( cube.nspec, 4*cube.n[2]+1 ) )
    #
    for iat in range( cube.nat ) :
        ispec = cube.spec_list.index( cube.atsymbol[iat] )
        for j in range( 4*cube.n[0]+1 ) : 
            spec_x[ispec,j] += gaussian( x_axis[j], cube.atpos[iat,0], 0.25) 
        for j in range( 4*cube.n[1]+1 ) : 
            spec_y[ispec,j] += gaussian( y_axis[j], cube.atpos[iat,1], 0.25) 
        for j in range( 4*cube.n[2]+1 ) : 
            spec_z[ispec,j] += gaussian( z_axis[j], cube.atpos[iat,2], 0.25) 
    #
    print("*** </ABUNDANCIES> ***")
    print(" ")
    #
    print("*** <PLOT> ***")
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
#
    fig = plt.figure(1)
    fig.subplots_adjust(hspace=.5)
# 
#   x
#
    ax1 = fig.add_subplot(3,1,1)
    for ispec in range(cube.nspec) : 
        plt.plot(x_axis,spec_x[ispec,:],label=atomic_symbol[int(cube.spec_list[ispec])])
    plt.xlim([np.min(x_axis[:]),np.max(x_axis[:])])
    plt.ylim([np.min(spec_x[:,:]),np.max(spec_x[:,:])])
    plt.xlabel('x (a.u.)')
    plt.ylabel('Abund.')
    plt.legend(loc=1,prop={'size':6})
#
# 
#   y
#
    ax2 = fig.add_subplot(3,1,2)
    for ispec in range(cube.nspec) : 
        plt.plot(y_axis,spec_y[ispec,:],label=atomic_symbol[int(cube.spec_list[ispec])])
    plt.xlim([np.min(y_axis[:]),np.max(y_axis[:])])
    plt.ylim([np.min(spec_y[:,:]),np.max(spec_y[:,:])])
    plt.xlabel('y (a.u.)')
    plt.ylabel('Abund.')
    plt.legend(loc=1,prop={'size':6})
# 
#   z
#
    ax3 = fig.add_subplot(3,1,3)
    for ispec in range(cube.nspec) : 
        plt.plot(z_axis,spec_z[ispec,:],label=atomic_symbol[int(cube.spec_list[ispec])])
    plt.xlim([np.min(z_axis[:]),np.max(z_axis[:])])
    plt.ylim([np.min(spec_z[:,:]),np.max(spec_z[:,:])])
    plt.xlabel('z (a.u.)')
    plt.ylabel('Abund.')
    plt.legend(loc=1,prop={'size':6})
# 
    plt.savefig(output_file_name)
    print("output written in : ", output_file_name)
    print("waiting for user to close image preview...")
    plt.show()
    fig.clear()
    print("*** </PLOT> ***")
    print(" ")
    #
    footer()

##
## Plot LDOS
##

def plot_ldos(energy_file_name,output_file_name,erange) : 
    #
    ### USER DEFINITION ###
    #energy_file_name      = 'energy.dat'
    #output_file_name      = 'ldos.png'
    #erange                = [40.,60.,0.01]
    ### END ##############
    #
    # N.B. 
    # the file energy_file_name should contain 4 values per line : en, sigma, weight, name
    # erange[0] = emin 
    # erange[1] = emax 
    # erange[2] = deltae 
    #
    import numpy as np
    import scipy as sp
    import matplotlib as mpl
    #
    header()
    versions()
    #
    print("*** <INPUT> ***")
    print("Energy file name     : ", energy_file_name)
    print("Output file name     : ", output_file_name)
    print("erange               : ", erange)
    print("*** </INPUT> ***")
    print(" ")
    #
    print("*** <LDOS> ***")
    #
    npte = int((erange[1]-erange[0])/erange[2]) + 1
    print("Using npte : ", npte) 
    # 
    en = EnergyFile(energy_file_name)
    en.show_content()
    #
    e_axis = np.linspace(erange[0], erange[1], npte, endpoint=True)
    #
    print("Requested (emin,emax) : ", erange[0], erange[1])
    print("Detected  (emin,emax) : ", np.min(en.e[:,0]), np.max(en.e[:,0]))
    #
    en.crunch_datain()
    #
    ntimes = 3
    #
    x_axis = np.linspace(0., np.max(en.axis[:]), ntimes*en.npt, endpoint=True)
    #
    xmap = np.zeros( ( ntimes*en.npt   , npte   ) )
    ymap = np.zeros( ( ntimes*en.npt   , npte   ) )
    cmap = np.zeros( ( ntimes*en.npt-1 , npte-1 ) )
    #
    for i in range(ntimes*en.npt) :
        for j in range(npte) :
            xmap[i,j] = x_axis[i]
            ymap[i,j] = e_axis[j]
    #
    from scipy import interpolate
    for ie in range(en.ne) :
        print(ie+1, "/", en.ne)
        for i in range(ntimes*en.npt-1) :
            tck = interpolate.splrep(en.axis[:], en.plot[ie,:], s=0, per=True)
            f   = interpolate.splev(x_axis, tck)
            for j in range(npte-1) :
                cmap[i,j] += f[i] * gaussian( e_axis[j], en.e[ie,0], en.e[ie,1] ) * en.e[ie,2] 
    #
    print("*** </LDOS> ***")
    print(" ")
    #
    print("*** <PLOT> ***")
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
#
    fig = plt.figure(1)
#
    ax = fig.add_subplot(1,1,1)
    plt.xlabel('avg. dir. (a.u.)')
    heatmap = ax.pcolormesh(xmap,ymap,cmap)
#    cbar = plt.colorbar(heatmap, orientation='horizontal', shrink=0.5, aspect=40)
    plt.xlim([0,np.max(x_axis[:])])
    plt.ylim([erange[0],erange[1]])
    plt.ylabel('Energy (eV)')
    plt.title("LDOS")
#    plt.axis('image')
#
#
    plt.savefig(output_file_name)
    print("output written in : ", output_file_name)
    print("waiting for user to close image preview...")
    plt.show()
    fig.clear()
    print("*** </PLOT> ***")
    print(" ")
    #
    footer()
