#!/usr/bin/python3

import cgi, cgitb
import sys, re
import math

form = cgi.FieldStorage()

print("Content-type: text/plain")
print()

###############################################################################
# User Input
###############################################################################

try:
 float(form.getvalue('B_MHZ'))
 B_MHZ = float(form.getvalue('B_MHZ'))
 B_MHZs = str(form.getvalue('B_MHZ'))
except ValueError:
 print("ERROR: Please enter a single value for the B0 field strength.")
 sys.exit()

try:
 float(form.getvalue('S2'))
 S2 = float(form.getvalue('S2'))
 S2s = str(form.getvalue('S2'))
except ValueError:
 print("ERROR: Please enter a single value for the order parameter.")
 sys.exit()

try:
 float(form.getvalue('tauR'))
 tauR = float(form.getvalue('tauR'))
 tauRs = str(form.getvalue('tauR'))
except ValueError:
 print("ERROR: Please enter a single value for the correlation time.")
 sys.exit()

try:
 float(form.getvalue('tauE'))
 tauE = float(form.getvalue('tauE'))
 tauEs = str(form.getvalue('tauE'))
except ValueError:
 print("ERROR: Please enter a single value for the correlation time.")
 sys.exit()

try:
 float(form.getvalue('distance'))
 distance = float(form.getvalue('distance'))
 distances = str(form.getvalue('distance'))
except ValueError:
 print("ERROR: Please enter a single value for the distance.")
 sys.exit()

try:
 float(form.getvalue('SE'))
 SE = float(form.getvalue('SE'))
 SEs = str(form.getvalue('SE'))
except ValueError:
 print("ERROR: Please enter a single value for the electron spin number.")
 sys.exit()

try:
# float(form.getvalue('nucleus_type'))
 form.getvalue('nucleus_type')
 nucleus_type = str(form.getvalue('nucleus_type'))
 nucleus_types = str(form.getvalue('nucleus_type'))
except ValueError:
 print("ERROR: Please enter a single value for the nucleus type.")
 sys.exit()

###############################################################################
# Log File
###############################################################################

from os import getenv
from datetime import datetime

now = datetime.now()

ipaddr = (getenv("HTTP_CLIENT_IP") or
getenv("HTTP_X_FORWARDED_FOR") or
getenv("REMOTE_ADDR") or
"UNKNOWN")

logfile = open('../logs/pre-calc.log', 'a')
logfile.write(now.strftime("%Y-%m-%d %H:%M"))
logfile.write('\t')
logfile.write(ipaddr)
logfile.write('\t')
logfile.write(B_MHZs)
logfile.write('\t')
logfile.write(distances)
logfile.write('\t')
logfile.write(SEs)
logfile.write('\t')
logfile.write(tauRs)
logfile.write('\t')
logfile.write(tauEs)
logfile.write('\t')
logfile.write(S2s)
logfile.write('\t')
logfile.write(nucleus_types)
logfile.write('\n')
logfile.close()

###############################################################################
# Constants
###############################################################################

# Pi
PI = math.pi

# mu-naught (permeability of free space, m kg s^-2 A*-2)
MU_0 = 1.25663706e-6

# mu-B (J T^-1)
MU_B = 9.2741e-24

# gyromagnetic ratio of proton (T^-1 s^-1)
GH = 42.5774821 * 2.0 * PI * 1.0e6
GN = -4.3156 * 2.0 * PI * 1.0e6
GC = 10.7084 * 2.0 * PI * 1.0e6

# electron g-factor (dimensionless)
GE = 2.0

###############################################################################
# Derived Constants
###############################################################################

# proton/nitrogen frequencies (rad/s)
WH = B_MHZ * 2.0 * PI * 1.0e6
WN = WH * GN / GH
WC = WH * GC / GH

# distance ^ -6
dist = distance * 1.0e-10
distpow = math.pow(dist,-6.0)

if nucleus_type == "1H":
    GX = GH
    WX = WH
if nucleus_type == "13C":
    GX = GC
    WX = WC
if nucleus_type == "15N":
    GX = GN
    WX = WN

###############################################################################
# Functions and  Observables
###############################################################################

def R1_PRE(J):
    return (2.0/5.0)*(MU_0/4.0/PI)*(MU_0/4.0/PI)*GX*GX*GE*GE*MU_B*MU_B*SE*(SE+1.0)*J(WX)

def R2_PRE(J):
    return (1.0/15.0)*(MU_0/4.0/PI)*(MU_0/4.0/PI)*GX*GX*GE*GE*MU_B*MU_B*SE*(SE+1.0)*(4.0*J(0.0)+3.0*J(WX))

###############################################################################
# Spectral Density Functions
###############################################################################

def J_mod1(tau):
    def J(w):
        return distpow*S2*tau/(1.0 + w*w*tau*tau)
    return J

tR = tauR + 0.0000001
tE = tauE + 0.0000001

tRp = tR * 1.0e-9
tEp = tE * 1.0e-9

tCp = 1.0/(1.0/tRp + 1.0/tEp)

Jp = J_mod1(tCp)

PRE_R1 = R1_PRE(Jp)
PRE_R2 = R2_PRE(Jp)

tauC = tCp * 1.0e9

###############################################################################
# Output
###############################################################################i

print("Script to calculate PRE rates")
print("(simple version)")
print()
print("Reference:")
print("Anthis N.J. & Clore G.M. (2014) Visualizing dark states by NMR spectroscopy,")
print("Quarterly Reviews of Biophysics, in press.")
print()
print("Input:")
print("field\t%s\tMHz" % B_MHZ)
print("dist.\t%s\tA" % distance)
print("e_spin\t%s\t" % SE)
print("tauC\t%.3f\tns" % tauC)
print("tauR\t%.3f\tns" % tauR)
print("tauE\t%.3f\tns" % tauE)
print("S2\t%s\t" % S2)
print("nucleus\t%s\t" % nucleus_type)
print()
print("Results:")
print("PRE_R1\t%13.5e\ts-1" % (PRE_R1))
print("PRE_R2\t%13.5e\ts-1" % (PRE_R2))

