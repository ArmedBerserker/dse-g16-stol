import numpy as np
#from Raymer p 517 or part 5 class 1
#variables
kgtolbs=2.20462 #conversion
slugtokg=1.35581795#slug*feet2 to kg*m2
mtof=3.28084 #meter to feet
b=14.81*mtof #span [feet]
g=32.2# [ft/s2]
L= 11*mtof#airplane length [feet]
r_cg=np.array([1,1,1])#[feet] from ground and axis of symmetry nose
m= 1809*kgtolbs# [pounds]
#approach
def Inertia_Class1():
    Rx=[0.260,0.251,0.373,0.240]
    Ry=[0.329,0.327,0.313,0.269]
    Rz=[0.399,0.391,0.384,0.461]
    Rx=np.mean(Rx)
    Ixx=(b**2*m*Rx**2)/(4*g)
    Ry=np.mean(Ry)
    Iyy=(L**2*m*Ry**2)/(4*g)
    Rz=np.mean(Rz)
    Izz=((b+L)/(2))**2*(m*Rz**2)/(4*g) #checked
    return(Ixx,Iyy,Izz)


#for rudder sizing we need it around the landing gear
def Inertia_parallel(I_cg,r,x_cg):
#all dimesions defined from ground and axis of symmetry nose
    d=x_cg-r
    x,y,z=d
    shift = m * np.array([
        [y**2+z**2,  -x*y,       -x*z],
        [-x*y,       x**2 + z**2,  -y*z],
        [-x*z,       -y*z,       x**2+y**2]
    ])
    return I_cg+shift



#class 2 part 5, p121
#symmetrical airplane thus Ixy and Iyz are 0
#self inertia of each component
#https://digitalcommons.usu.edu/cgi/viewcontent.cgi?article=1030&context=mae_stures


def self_inertia_cylinder(m, length, radius):
    #engine
    Ixx_i = 0.5 * m/g* radius ** 2
    Iyy_i = m/g * (3 * radius ** 2 + length ** 2) / 12
    Izz_i = Iyy_i
    return (Ixx_i, Iyy_i, Izz_i)

def self_inertia_fuselage(m, length, radiusinterior,radiusexterior):
    #fuselage use hollow cylinder for better approximation

    Ixx_i = 0.5 * m/g * (radiusexterior ** 2+radiusinterior**2)
    Iyy_i = m/g* (3 *(radiusexterior ** 2+radiusinterior**2)+length**2 ) / 12
    Izz_i = Iyy_i
    return (Ixx_i, Iyy_i, Izz_i)


def self_inertia_full_rectangle(m, lx, ly, lz):
    #cargo, capacitor, fuel tanks, wing, stabiliser,...
    #check that the wing is not hollow

    Ixx_i = m/g * (ly ** 2 + lz ** 2) / 12
    Iyy_i = m/g * (lx ** 2 + lz ** 2) / 12
    Izz_i = m/g * (lx ** 2 + ly ** 2) / 12
    return (Ixx_i, Iyy_i, Izz_i)


#inertia is 0 for lumped mass such as fuel system, electronics,... assumed to be concentrated mass somewhere


def Inertia_Class2(components, xcg, ycg, zcg):
    Ixx = Iyy = Izz = 0.0
    Ixy = Iyz = Ixz = 0.0

    for c in components:
        #find distance for each component
        m_i = c['m_i']
        dx = c['x_i']*mtof- xcg
        dy = c['y_i']*mtof- ycg
        dz = c['z_i']*mtof - zcg

        #add to inertia parallel theorem
        Ixx += m_i/g * (dy ** 2 + dz ** 2)
        Iyy += m_i/g * (dz ** 2 + dx ** 2)
        Izz += m_i/g * (dx ** 2 + dy ** 2)
        Ixy += m_i/g * (dx * dy)
        Iyz += m_i/g * (dy * dz)
        Ixz += m_i/g * (dz * dx)

        if c.get('name') == 'Wing':
            m_i=c['m_i']
            lx=0.85*mtof
            ly=7.1868*mtof
            lz=0.102*mtof
            Ixx_rec, Iyy_rec, Izz_rec = self_inertia_full_rectangle(m_i, lx, ly, lz)#change when have number

            Ixx += Ixx_rec
            Iyy += Iyy_rec
            Izz += Izz_rec
        if c.get('name') == 'Fuselage':
            m_i = c['m_i']
            length=11*mtof
            radiusinterior=0.72*mtof
            radiusexterior=0.725*mtof
            Ixx_rec, Iyy_rec, Izz_rec = self_inertia_fuselage(m_i, length, radiusinterior,radiusexterior)

            Ixx += Ixx_rec
            Iyy += Iyy_rec
            Izz += Izz_rec
        if c.get('name') == 'Verticaltail':
            m_i = c['m_i']
            lx = 0.9*mtof
            ly =0.126*mtof
            lz =1.2*mtof
            Ixx_rec, Iyy_rec, Izz_rec = self_inertia_full_rectangle(m_i, lx, ly, lz)#change when have number

            Ixx += Ixx_rec
            Iyy += Iyy_rec
            Izz += Izz_rec

        if c.get('name') == 'Horizontaltail':
            m_i = c['m_i']
            lx = 0.4*mtof
            ly = 2.26*mtof
            lz = 0.05*mtof
            Ixx_rec, Iyy_rec, Izz_rec = self_inertia_full_rectangle(m_i, lx, ly, lz)#change when have number

            Ixx += Ixx_rec
            Iyy += Iyy_rec
            Izz += Izz_rec
        if c.get('name') == 'Engine':
            m_i=c['m_i']
            length=1.1*mtof
            radius=0.315*mtof
            Ixx_rec, Iyy_rec, Izz_rec = self_inertia_cylinder(m_i, length, radius)

            Ixx += Ixx_rec
            Iyy += Iyy_rec
            Izz += Izz_rec

    return dict(Ixx=Ixx, Iyy=Iyy, Izz=Izz,
                Ixy=Ixy, Iyz=Iyz, Ixz=Ixz)


#call everything
#inertia class 1

#add all the components given in class 2
components = [
    # name           m_i [lb]   x_i    y_i    z_i   [meter] shape
    {'name': 'Wing', 'm_i': 155.9*kgtolbs, 'x_i': 6.212, 'y_i': 0.0, 'z_i': 2.96}, #good mass
    {'name': 'Fuselage', 'm_i': 174.2*kgtolbs, 'x_i': 5.5, 'y_i':0.0 , 'z_i': 2.11},#good mass
    {'name': 'Verticaltail', 'm_i': 17*kgtolbs, 'x_i': 8.01, 'y_i': 0.0, 'z_i': 3.76}, #good mass
    {'name': 'Horizontaltail', 'm_i': 27.6*kgtolbs, 'x_i': 8.01, 'y_i': 0.0, 'z_i':2.96}, #good mass
    {'name': 'Engine', 'm_i': 82.2*2*kgtolbs, 'x_i': 6.04, 'y_i': 0.0, 'z_i': 2.435}, #good mass +nacelle mass
    {'name': 'LandingGear1', 'm_i': 9*kgtolbs, 'x_i': 0.605, 'y_i': 0.0, 'z_i': 0.63},  #good
    {'name': 'LandingGear2', 'm_i': 15*2*kgtolbs, 'x_i': 6.52, 'y_i': 0.0, 'z_i': 0.63},  #good
    {'name': 'Fuel', 'm_i': (167.8+54)*kgtolbs, 'x_i':6.212, 'y_i': 0.0, 'z_i': 2.96}, #good
    {'name': 'Cargo', 'm_i': 200 * kgtolbs, 'x_i': 6.7+0.25, 'y_i': 0.0, 'z_i': 2.06},#good
    {'name': 'Row1', 'm_i': 154 * kgtolbs, 'x_i': 2.92+0.25, 'y_i': 0.0, 'z_i': 1.91},#good
    {'name': 'Row2', 'm_i': 154 * kgtolbs, 'x_i': 5+0.25, 'y_i': 0.0, 'z_i': 1.91},#good
    {'name': 'Row3', 'm_i': 154 * kgtolbs, 'x_i': 6+0.25, 'y_i': 0.0, 'z_i': 1.91},#good
    {'name': 'Supercapacitor', 'm_i': 154 * kgtolbs, 'x_i': 6 + 0.25, 'y_i': 0.0, 'z_i': 1.91},  # good


]


total_W = sum(c['m_i'] for c in components)
xcg_calc = sum(c['m_i'] * c['x_i']*mtof for c in components) / total_W
ycg_calc = sum(c['m_i'] * c['y_i']*mtof  for c in components) / total_W
zcg_calc = sum(c['m_i'] * c['z_i']*mtof  for c in components) / total_W

print(f"\nClass II  computed c.g.:  x={xcg_calc:.3f} ft,  "
      f"y={ycg_calc:.3f} ft,  z={zcg_calc:.3f} ft")


Ixx, Iyy, Izz = Inertia_Class1()
print(f"Inertia class 1 Ixx= {Ixx:.3f},Iyy={Iyy:.3f},Izz={Izz:.3f}")
print(f"Inertia class 1 normal Ixx= {Ixx*slugtokg:.3f},Iyy={Iyy*slugtokg:.3f},Izz={Izz*slugtokg:.3f}")
#parallel set  up
x_cg=[xcg_calc,ycg_calc,zcg_calc]
r_lg = np.array([1, 2, 3]) #example for landing gear [feet]
I_cg=np.array([
    [Ixx, 0, 0],
    [0, Iyy, 0],
    [0, 0, Izz]
])



I=Inertia_parallel(I_cg,r_lg,x_cg)
I2 = Inertia_Class2(components, xcg_calc, ycg_calc, zcg_calc)

print(f"\nClass II  Ixx = {I2['Ixx']:.3f} slug·ft²")
print(f"             Iyy = {I2['Iyy']:.3f} slug·ft²")
print(f"             Izz = {I2['Izz']:.3f} slug·ft²")
print(f"             Ixy = {I2['Ixy']:.3f} slug·ft² ")
print(f"             Iyz = {I2['Iyz']:.3f} slug·ft² ")
print(f"             Ixz = {I2['Ixz']:.3f} slug·ft²")

# Full inertia tensor
I2 = np.array([
    [I2['Ixx'], -I2['Ixy'], -I2['Ixz']],
    [-I2['Ixy'], I2['Iyy'], -I2['Iyz']],
    [-I2['Ixz'], -I2['Iyz'], I2['Izz']],
])
print("\nClass II inertia tensor [slug·ft²]:")
print(np.round(I2, 3))
#check if Ixy and Iyz are close to 0