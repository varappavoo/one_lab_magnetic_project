#!/usr/bin/env python
#coding=utf-8

# 1) Launch "gmsh pend.py", or open "pend.py" with Gmsh's File->Open menu
# 2) Click on "Run" in the left Gmsh panel
# 3) there is no number 3... :-)

# To interact with ONELAB, the script must first import the onelab.py module:
import onelab
import math, os

# The script then creates a ONELAB client with:
c = onelab.client(__file__)

# Creating the client connects the script to the onelab server, through a
# socket. The __file__ argument is a python variable. It tells ONELAB in which
# directory the script being executed is located.

def exportMsh(le1,le2):
   mshFile = open(c.getPath("pend.msh"), 'w')
   mshFile.write('$MeshFormat\n2.2 0 8\n$EndMeshFormat\n')
   mshFile.write('$Nodes\n3\n1 0 0 0\n2 0 %s 0\n3 0 %s 0\n$EndNodes\n' %(-le1, -le1-le2))
   mshFile.write('$Elements\n3\n1 1 2 0 1 1 2\n2 1 2 0 1 2 3\n3 15 2 0 2 3\n$EndElements\n')
   mshFile.close()

def exportMshOpt():
   optFile = open(c.getPath("pend.msh.opt"),'w')
   optFile.write('n = PostProcessing.NbViews - 1;\n')
   optFile.write('If(n >= 0)\nView[n].ShowScale = 0;\nView[n].VectorType = 5;\n')
   optFile.write('View[n].ExternalView = 0;\nView[n].DisplacementFactor = 1 ;\n')
   optFile.write('View[n].PointType = 1;\nView[n].PointSize = 5;\n')
   optFile.write('View[n].LineWidth = 2;\nEndIf\n')
   optFile.close()

def exportIter(iter,t,x1,y1,x2,y2):
   mshFile = open(c.getPath("pend.msh"),'a')
   mshFile.write('$NodeData\n1\n"motion"\n1\n\t%f\n3\n\t%d\n3\n' % (t, iter))
   mshFile.write('\t3\n\t1 0 0 0\n\t2 %f %f 0\n\t3 %f %f 0\n$EndNodeData\n' %(x1,y1,x2,y2))
   mshFile.close()

g = 9.8	# acceleration of gravity
m = 0.3 # mass of pendulum balls

# New ONELAB variables can then be defined using defineNumber, e.g.:
l = c.defineNumber('Geom/arm length [m]', value=1.0)
time = c.defineNumber('Dyna/time [s]', value=0.0)
dt = c.defineNumber('Dyna/time step [s]', value=0.001)
tmax = c.defineNumber('Dyna/max time [s]', value=20)
refresh = c.defineNumber('Dyna/refresh interval [s]', value=0.1)
theta0 = c.defineNumber('Init/initial theta angle [deg]', value=10, 
                         attributes={'Highlight':'Pink'})
phi0 = c.defineNumber('Init/initial phi angle [deg]', value=180,
                       attributes={'Highlight':'Pink'})

# When the script is run, if the parameter Geom/arm length [m] has not been
# previously defined, it takes the value (1.0) provided in defineNumber and is
# sent to the ONELAB server. The "/" character in the variable name is
# interpreted as a path separator, and results in the creation of a sub-tree in
# the graphical user interface. If the script is re-run later, the value will be
# updated using the value from the server (unless it is labeled as readOnly: see
# below). When Gmsh runs a ONELAB client, the client can be run in two modes:
# c.action=='check' to check the coherence of the ONELAB database and make
# adjustments if necessary, and c.action=='compute' to perform the actual
# computation. For instance, in 'check' mode, the double pendulum client simply
# defines the ONELAB variables it wants to share with the server, then exits:
if c.action == 'check' :
   exit(0)

# In 'compute' mode, the code enters a loop and performs the actual computation.
   
l1 = l;
l2 = l;
m1 = m;
m2 = m;
theta = theta0 / 180.*math.pi;
phi = phi0 / 180.*math.pi;
theta_dot = 0.0
phi_dot = 0.0
refr = 0.0
iter = 0
time = 0.0

while (time < tmax):
   delta = phi - theta
   sdelta = math.sin(delta)
   cdelta = math.cos(delta)
   theta_dot_dot = ( m2*l1*(theta_dot**2.0)*sdelta*cdelta
                     + m2*g*math.sin(phi)*cdelta
                     + m2*l2*(phi_dot**2.0)*sdelta
                     - (m1+m2)*g*math.sin(theta) )
   theta_dot_dot /= ( (m1+m2)*l1 - m2*l1*(cdelta)**2.0 )
   
   phi_dot_dot = ( -m2*l2*(phi_dot**2.0)*sdelta*cdelta
                    + (m1+m2)*(g*math.sin(theta)*cdelta
                               - l1*(theta_dot**2.0)*sdelta
                               - g*math.sin(phi)) )
   phi_dot_dot /= ( (m1+m2)*l2 - m2*l2*(cdelta)**2.0 )
   
   theta_dot = theta_dot + theta_dot_dot*dt
   phi_dot = phi_dot + phi_dot_dot*dt

   theta = theta + theta_dot*dt
   phi = phi + phi_dot*dt

   x1 =  l1*math.sin(theta)
   y1 = -l1*math.cos(theta)
   x2 =  l1*math.sin(theta) + l2*math.sin(phi)
   y2 = -l1*math.cos(theta) - l2*math.cos(phi)

   time += dt
   refr += dt

   exportMshOpt()

   if refr >= refresh:
      refr = 0
      # During the computation the script can directly set a value in the ONELAB
      # database with setNumber:
      c.setNumber(c.name + '/Progress', value=time, min=0, max=tmax, visible=0)
      c.setNumber('Dyna/time [s]', value=time)
      c.setNumber('Solu/phi', value=phi)
      c.addNumberChoice('Solu/phi', phi)
      c.setNumber('Solu/theta', value=theta)
      c.addNumberChoice('Solu/theta', theta)
      c.setNumber('Solu/phi dot', value=phi_dot)
      c.addNumberChoice('Solu/phi dot', phi_dot)
      c.setNumber('Solu/theta dot', value=theta_dot)
      c.addNumberChoice('Solu/theta dot', theta_dot)

      # It can also ask the server to refresh...
      c.setString('Gmsh/Action', value='refresh')

      # or to stop...
      if(c.getString(c.name + '/Action') == 'stop'):
         break;

      exportMsh(l1, l2)
      exportIter(iter, time, x1, y1+l1, x2, y2+l1+l2)

      # or to read a file:
      c.mergeFile(c.checkPath('pend.msh'))
      iter += 1

      # The check path function (c.checkPath) builds the pathname of a file
      # named pend.msh located in the same directory as the script under
      # execution, and then checks on whether the pathname exists on disk. If
      # not, an error message is issued. Use the regular path function
      # (c.getPath) to build a pathname without checking on the presence on disk
      # of the file.

c.setNumber(c.name + '/Progress', value=0)
