from classes.aircraft_2 import loader, Aircraft
from class1 import prelim_drag
from class1.c1_m import Class1

aircraft = loader.load('yamls/aircraft.yaml', Aircraft)

thing = prelim_drag.prelim_drag(aircraft)
print(thing)

thing2 = Class1(aircraft)

print(thing2.fuel_frac_used())