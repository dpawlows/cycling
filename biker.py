# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 11:08:01 2023

@author: altha
"""

import numpy as np
import matplotlib.pyplot as pp
from scipy.signal import savgol_filter
from sys import exit
import gpxpy
import gpxpy.gpx
import geopy.distance

import numpy as np
filename = 'Alpe_d_Huez.gpx' 
gpx_file = open(filename)

gpx = gpxpy.parse(gpx_file) 
h = []
latitude = []
longitude = []
for track in gpx.tracks:
        for segment in track.segments:
            for point in segment.points:
                print:'Point at ({0},{1}) -> {2}'.format(point.latitude, point.longitude, point.elevation)
                h.append(point.elevation)
                latitude.append(point.latitude)
                longitude.append(point.longitude)
                if len(latitude) > 2:
                    if latitude[-1] == latitude[-2]:
                        #if this happens, we found a duplicate point in the gpx file
                        #so just remove it from the lists
                        latitude.pop()
                        longitude.pop()
                        h.pop()
                     
h = np.array(h)
hSmoothed = savgol_filter(h, 200, 3) #Smooth the elevation data a little

elevationEarth = h - h[0]

npoints = len(h)
latitude = np.array(latitude)
dlat = latitude[1:] - latitude[0]

longitude = np.array(longitude)
dlong = longitude[1: ] - longitude[0]


distanceTraveled = [0.0]
angle = []
previousCoords = (latitude[0],longitude[0])
gradient = []
for i in range(1,npoints):
    #Set the distance variable as the distance traveled along the GPS
    #We will use this array for interpolation when solving for the motion
    #of the biker.

    coords = (latitude[i],longitude[i]) #new coordinates
    #how far between current and previous coordinates
    newMovement = geopy.distance.geodesic(previousCoords,coords).m 

    #keep track
    distanceTraveled.append(distanceTraveled[-1] + newMovement)

    #get change in elevation
    y = h[i] - h[i-1]

    #gradient is rise over run
    gradient.append(y/newMovement)

    #angle is arctan of gradient
    angle.append(np.arctan(gradient[-1]))

    previousCoords = coords #new coords become the old coords for the next iteration    


# Set constants________________________________________________________________
v01 = 4 #initial velocity   
def P1():
    p1 = 400
    return p1 #generate a random power output for dt 
m = 70 #mass of biker (Kg)

dt = 0.1 #time step
final_time = 200 #total time 
step_value = int(final_time / dt) #amount of iterations


C_d = 0.9 #coefficient of drag
A = 0.33 # Area of biker (m^2)
ro = 1.225 #density of air (Kg/m^3)
n = 2E-5 # viscosity of air (Pa*s)
heightBiker = 1.6 #m


g = 9.81 #accerlation due to gravity (m/s^2)
C_r = 0.008 #coefficient of rolling resitance for a bike
C_h = 0.015 #coefficient of kinetic friction do to the wheel hub 
# set velocity_________________________________________________________________
v1 = [v01]
time = [0] #Need to track time
bikerElevation = [h[0]] #track elevation of the biker


x1 = [0]

#Don't go past the gps track!
endOfRoute = distanceTraveled[-1]

#Run the variables for every iteration________________________________________
while x1[-1] < endOfRoute:
    #First, get the angle of the route by interpolating the 
    #gps track data to the current position of the rider
    current_angle = np.interp(x1[-1],distanceTraveled[1:],angle)
    bikerElevation.append(np.interp(x1[-1],distanceTraveled,h))

    Power_output = P1() #update the power

    #calculate the forces (units are acceleration!!!!)
    cyclist = Power_output/(m*v1[-1]) 
    gravitationalForce = -g * np.sin(current_angle)
    normal_Force = g*np.cos(current_angle)
    Rolling_resistance = -C_r * normal_Force/m
    Air_drag = -0.5*C_d*ro*A*v1[-1]**2/m
    Viscous_force = -n*(v1[-1]/heightBiker*m)
    Internal_friction = -C_h * normal_Force/m

    #F = ma
    dvdt = (cyclist + Air_drag + Viscous_force + gravitationalForce + Rolling_resistance + Internal_friction) / (m)

    #Take the euler step to update the velocity and position of the biker
    v1.append(v1[-1] + dvdt*dt)
    x1.append(x1[-1] + v1[-1]*dt)
    #keep track of time
    time.append(time[-1]+dt)
    
# Plot Everything______________________________________________________________

vKPH = np.asarray(v1) * 3600/1000 #KPH are appropriate units for cycling
pos = filename.rfind('.gpx')
label = filename[:pos]
pp.figure(figsize=(8,12))
pp.subplot(311)
pp.plot(time, bikerElevation, label=label, color='b')
# pp.xlabel('Time (s)')
pp.ylabel('Elevation (m)')
pp.legend(frameon=False)


pp.subplot(312)
pp.plot(time, vKPH, color='b')
# pp.xlabel('Time (s)')
pp.ylabel('Velocity (kph)')

pp.subplot(313)
pp.plot(time, np.asarray(x1)/1000., label='Total Time: {} min'.format(int(time[-1]/60.)), color='b')
pp.xlabel('Time (s)')
pp.ylabel('Distance (km)')
pp.legend(frameon=False)

pp.savefig('plot.png')

