import numpy as np
from numpy import sqrt, sin, cos, pi
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from subprocess import Popen, PIPE

#########################################
##### Installs/Runs the C++ program #####
#########################################

def compile_cpp_files(directory, filenames, compiler = "g++", flags = "-o"):
    for filename in filenames:
        print( "Compiling", filename + ".cpp" )
        p = Popen( [ compiler, flags, directory+filename, directory+filename+".cpp"] , stdout = PIPE)
        out = p.communicate()
        print( "Output:", out[0], "\nErrors:", str(out[1]), "\n" )


def run_treecode(path, directory, inputfile, outputfiles, dt, nsteps, eps, threshold, xmin=-100, xmax=100, ymin=-100, ymax=100, zmin=-100, zmax=100):
    p = Popen( [path, directory, inputfile, outputfiles, str(dt), str(nsteps), str(eps), str(threshold), str(xmin), str(xmax), str(ymin), str(ymax), str(zmin), str(zmax) ] , stdout = PIPE)
    out = p.communicate()
    #print( "Output:", out[0], "\nErrors:", str(out[1]), "\n" )
    return 


########################################
##### Data management and plotting #####
########################################

def combine_galaxies(directory, target, companion, offset = [ 0, 0, 10.0 ], vel = [ 0, 0, 0 ], angle = [ 0., 0.7] ):
    filename = target + "_" + companion + ".txt"
    with open(directory+"data/" + target, "r") as f:
        lines = []    
        for i, line in enumerate(f.read().splitlines()):
            if i == 0:
                lines.append( int(line.split()[0]) )
            else:
                numbers = line[1:].split()
                numbers = [ float(n) for n in numbers ]
                numbers = [ str(n) for n in numbers ]
                lines.append( numbers[0] )
                for n in numbers[1:]:
                    lines[-1] += "," + n
    if len(companion) > 1: 
        with open(directory+"data/" + companion, "r") as f:
            for i, line in enumerate(f.read().splitlines()):
                if i == 0:
                    lines[0] += int( line.split()[0] )
                else:
                    numbers = line[1:].split()
                    numbers = [ float(n) for n in numbers ]
                    for j in range(3):
                        if offset[j] != 0: 
                            numbers[1+j] += offset[j]
                        if vel[j] != 0:
                            numbers[4+j] += vel[j]
                    numbers = [ str(n) for n in numbers ]
                    lines.append( numbers[0] )
                    for n in numbers[1:]:
                        lines[-1] += "," + n
    lines[0] = "Parameters:" + str(lines[0])
    with open(directory+"data/" + filename, "w+") as f:
        for i, line in enumerate(lines):
            if i != 0: f.write("\n")
            f.write(line)
    return "data/" + filename
    
def plot_3D(source, destination, groups, filerange = range(1,2000, 2), angle = [90., 0.], colors = [ 'lightyellow', 'violet', 'lightyellow', 'cyan', 'red' ], xlim=5, ylim=5, zlim = 5, axes = True):
    fig = plt.figure()
    fig.set_size_inches(10,10)
    ax = fig.gca(projection='3d')
    padding = len(str( len(filerange)-1) )
    for j, n in enumerate(filerange):
        with open(source + "_" + str(n) + ".txt", "r") as f:
            time = -1
            x, y, z = [], [], []
            for i, line in enumerate(f.read().splitlines()):
                if i==0:
                    time = float( line.split("=")[1] )
                else:
                    linedata = line.split(",")
                    x.append( float(linedata[0]) )
                    y.append( float(linedata[1]) )
                    z.append( float(linedata[2]) )
            ax.cla()
            a, b = 0, 0
            for i in range(len(groups)):
                b += groups[i]
                ax.scatter(x[a:b], y[a:b], z[a:b], c=colors[i], s=1, alpha=0.3) 
                a += groups[i]
            ax.set_aspect("equal")
            ax.set_xlim(-xlim,xlim)
            ax.set_ylim(-ylim,ylim)
            ax.set_zlim(-zlim,zlim)
            if axes:
                ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
                ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
                ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
                ax.w_xaxis.line.set_color("white")
                ax.tick_params(colors='white')
                ax.set_xlabel("x Distance (3 kpc)", color='white')
                ax.set_ylabel("y Distance (3 kpc)", color='white')
                ax.set_zlabel("z Distance (3 kpc)", color='white')
            else:
                ax.set_axis_off() 
            ax.set_facecolor('black')
            ax.text2D(0.05, 0.95, "Time: " + str(time) + " (" + str(time*10)[:6] + " million years)", transform=ax.transAxes, fontsize=15, color='white')
            if angle:
                ax.view_init(angle[0], angle[1])
            plt.savefig(destination + str(j).zfill(padding) + ".png" ) 


def plot_3D_fast(source, destination, groups, filerange = range(1,2000, 2), angle = [90., 0.], colors = [ 'lightyellow', 'violet', 'lightyellow', 'cyan', 'red' ], xlim=5, ylim=5, zlim = 5, axes = True):
    fig = plt.figure()
    fig.set_size_inches(10,10)
    ax = fig.gca(projection='3d')
    padding = len(str( len(filerange)-1) )
    for j, n in enumerate(filerange):
        with open(source + "_" + str(n) + ".txt", "r") as f:
            time = -1
            x, y, z = [], [], []
            for i, line in enumerate(f.read().splitlines()):
                if i==0:
                    time = float( line.split("=")[1] )
                else:
                    linedata = line.split(",")
                    x.append( float(linedata[0]) )
                    y.append( float(linedata[1]) )
                    z.append( float(linedata[2]) )
            ax.cla()
            a, b = 0, 0
            for i in range(len(groups)):
                b += groups[i]
                ax.scatter(x[a:b], y[a:b], z[a:b], c=colors[i], s=1, alpha=0.3, mode='2dcross') 
                a += groups[i]
            ax.set_aspect("equal")
            ax.set_xlim(-xlim,xlim)
            ax.set_ylim(-ylim,ylim)
            ax.set_zlim(-zlim,zlim)
            ax.set_axis_off() 
            ax.set_facecolor('black')
            ax.text2D(0.05, 0.95, "Time: " + str(time) + " (" + str(time*10)[:6] + " million years)", transform=ax.transAxes, fontsize=15, color='white')
            if angle:
                ax.view_init(angle[0], angle[1])
            plt.savefig(destination + str(j).zfill(padding) + ".png" ) 
