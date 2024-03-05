# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 11:08:01 2023

@author: altha
"""

import numpy as np
import matplotlib.pyplot as pp
from sys import exit
import gpxpy
import gpxpy.gpx
import geopy.distance

import numpy as np
gpx_file = open('Alpe_d_Huez.gpx')

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
                
h = np.array(h)
elevationEarth = h - h[0]

npoints = len(h)
latitude = np.array(latitude)
dlat = latitude[1:] - latitude[0]

longitude = np.array(longitude)
dlong = longitude[1: ] - longitude[0]

#position = (np.sin(dlat/2))**2 + np.cos(latitude) * np.cos(latitude) * (np.sin(dlong/2))**2
#c = 2 * np.atan2(np.sqrt(position), np.sqrt(1-position))
#distance = position * c


coords0 = (latitude[0],longitude[0])
distanceEarth = []
angle = []

for i in range(npoints):
    coords = (longitude[i],latitude[i])

    distanceEarth.append(geopy.distance.geodesic(coords0,coords).km)
    
for i in range(1,npoints):
    x = distanceEarth[i] - distanceEarth[i-1]
    y = elevationEarth[i] - elevationEarth[i-1]
    angle.append(np.arctan(y/x))
    
    
breakpoint()

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



g = 9.81 #accerlation due to gravity (m/s^2)
C_r = 0.008 #coefficient of rolling resitance for a bike
C_h = 0.015 #coefficient of kinetic friction do to the wheel hub 
# set velocity_________________________________________________________________
V1 = [v01]
time = np.linspace(0, final_time, step_value)

#______________________________________________________________________________

# Define functions for slopes that vary with time
#def Park_ride(t):
   # path = np.random.uniform(-5,5)
    #return np.sin(np.radians(path))
    

#def Hiking_trip(t):
    #return .1*np.sin(.01*t)  

#def X_games(t):
    #return .5*np.sin(.005*t) 

#print("1. Park Ride")
#print("2. Hiking Trip")
#print("3. X Games")

#try:
 #   routeType = int(input("Choose 1,2 or 3:      "))
#except:
 #   print("Invalid entry")
  #  exit(1)

#heights_Park_ride = [0]
#heights_Hiking_trip = [0]
#heights_X_games = [0]


#Run the variables for every itteration________________________________________
for i in range(1, step_value):
        
    Power_output = P1()

    normal_Force = m * g * np.sin(current_angle)
    
    Rolling_resistance = C_r * normal_Force
    
    v1 = V1[i - 1]
    
    Air_drag = 0.5*C_d*ro*A*v1**2
    
    Viscous_force = n*(v1/h)
    
    Internal_friction = C_h * normal_Force
    
    D1 = np.cumsum(np.array(V1)*dt)


    dvdt = (Power_output - Air_drag - Viscous_force - normal_Force - Rolling_resistance - Internal_friction) / (m)

    v_new = v1 + dvdt*dt
    x_new= D1 + (v_new*dt)
    
    current_angle = np.interp(x_new,distanceEarth,angle)
    
    V1.append(v_new)
    D1_total = np.sum(np.array(V1)*dt)
   

   

    print(D1_total)
    #print("Rolling force:",Rolling_resistance, "slope force:",slope_force, "power output:",Power_output, "drag force:",Air_drag)


average_velocity = np.mean(V1)
print("average velocity:",average_velocity)




# Plot Everything______________________________________________________________

#YOU IF YOU PLAN ON CHOSING A CERTAIN BIKE PATH YOU HAVE TO MODIFY LINE NUMBER 119 TO MATCH YOUR CHOISE BEFOREHAND
pp.plot(time, x_new, label='Dr.Pawlowski', color='orange')
pp.title('distance over Time')
pp.xlabel('Time (s)')
pp.ylabel('distance')
pp.legend()
pp.show()
