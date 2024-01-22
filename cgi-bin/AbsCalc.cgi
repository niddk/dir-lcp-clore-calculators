#!/usr/bin/python3

import cgi, cgitb
import sys, re

form = cgi.FieldStorage()

###############################################################################
# Protein Sequence
###############################################################################

seq1 = form.getvalue('sequence')

seq2 = re.sub('[^ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy]','',seq1)
seq = seq2.upper()

###############################################################################
# Log File
###############################################################################

from os import getenv
from datetime import datetime

now = datetime.now()

ipaddr = (getenv("HTTPS_CLIENT_IP") or
getenv("HTTPS_X_FORWARDED_FOR") or
getenv("REMOTE_ADDR") or
"UNKNOWN")

#logfile = open('/gmclore.org/public_html/weblogs/A205.log', 'a')
logfile = open('../logs/A205.log', 'a')
logfile.write(now.strftime('%Y-%m-%d %H:%M'))
logfile.write('\t')
logfile.write(ipaddr)
logfile.write('\t')
logfile.write(seq)
logfile.write('\n')
logfile.close()

###############################################################################
# Constants
###############################################################################

# Atomic masses
mwH = 1.007825
mwC = 12.0107
mwN = 14.0067
mwO = 15.9994
mwS = 32.065

# Side chain element counts
# non-exchangeable hydrogens
AHno = 3
CyHno = 2
DHno = 2
EHno = 4
FHno = 7
GHno = 1
HHno = 4
IHno = 9
KHno = 8
LHno = 9
MHno = 7
NHno = 2
PHno = 6
QHno = 4
RHno = 6
SHno = 2
THno = 4
VHno = 7
WHno = 7
YHno = 6
# exchangeable hydrogens
AHex = 0
CyHex = 1
DHex = 1 
EHex = 1
FHex = 0
GHex = 0
HHex = 1
IHex = 0
KHex = 2
LHex = 0
MHex = 0
NHex = 2
PHex = -1
QHex = 2
RHex = 4
SHex = 1
THex = 1
VHex = 0
WHex = 1
YHex = 1
# carbons
AC = 1
CyC = 1
DC = 2
EC = 3
FC = 7
GC = 0
HC = 4
IC = 4
KC = 4
LC = 4
MC = 3
NC = 2
PC = 3
QC = 3
RC = 4
SC = 1
TC = 2
VC = 3
WC = 9
YC = 7
# nitrogens
AN = 0
CyN = 0
DN = 0
EN = 0
FN = 0
GN = 0
HN = 2
IN = 0
KN = 1
LN = 0
MN = 0
NN = 1
PN = 0
QN = 1
RN = 3
SN = 0
TN = 0
VN = 0
WN = 1
YN = 0
# oxygens
AO = 0
CyO = 0
DO = 2
EO = 2
FO = 0
GO = 0
HO = 0
IO = 0
KO = 0
LO = 0
MO = 0
NO = 1
PO = 0
QO = 1
RO = 0
SO = 1
TO = 1
VO = 0
WO = 0
YO = 1
# sulfurs
AS = 0
CyS = 1
DS = 0
ES = 0
FS = 0
GS = 0
HS = 0
IS = 0
KS = 0
LS = 0
MS = 1
NS = 0
PS = 0
QS = 0
RS = 0
SS = 0
TS = 0
VS = 0
WS = 0
YS = 0

# Extinction coefficients at 280 nm in M-1 cm-1.
W280 = 5500
Y280 = 1490
CyCy280 = 125 # disulfide bond

# Extinction coefficient for each peptide bond at 205 nm in M-1 cm-1. 
BB205 = 2780

# Extinction coefficients for sidechains at 205 nm in M-1 cm-1.
W205 = 20400
F205 = 8600
Y205 = 6080
H205 = 5200
M205 = 1830
R205 = 1350
N205 = 400
Q205 = 400	
Cy205 = 690 
CyCy205 = 820 # disulfide bond: 2200 - 2*690 = 820.

###############################################################################
# Calculations
###############################################################################

A = seq.count('A') 
Cy = seq.count('C') 
D = seq.count('D') 
E = seq.count('E')
F = seq.count('F')
G = seq.count('G')
H = seq.count('H')
I = seq.count('I')
K = seq.count('K')
L = seq.count('L')
M = seq.count('M')
N = seq.count('N')
P = seq.count('P')
Q = seq.count('Q')
R = seq.count('R')
S = seq.count('S')
T = seq.count('T')
V = seq.count('V')
W = seq.count('W')
Y = seq.count('Y')

rescount = A + Cy + D + E + F + G + H + I + K + L + M + N + P + Q + R + S + T + V + W + Y
BB = rescount - 1
CyCy = Cy/2 

Ccount = 2*rescount + A*AC + Cy*CyC + D*DC + E*EC + F*FC + G*GC + H*HC + I*IC + K*KC + L*LC + M*MC + N*NC + P*PC + Q*QC + R*RC + S*SC + T*TC + V*VC + W*WC + Y*YC
Hexcount = 2 + 1*rescount + A*AHex + Cy*CyHex + D*DHex + E*EHex + F*FHex + G*GHex + H*HHex + I*IHex + K*KHex + L*LHex + M*MHex + N*NHex + P*PHex + Q*QHex + R*RHex + S*SHex + T*THex + V*VHex + W*WHex + Y*YHex
Hnocount = 1*rescount + A*AHno + Cy*CyHno + D*DHno + E*EHno + F*FHno + G*GHno + H*HHno + I*IHno + K*KHno + L*LHno + M*MHno + N*NHno + P*PHno + Q*QHno + R*RHno + S*SHno + T*THno + V*VHno + W*WHno + Y*YHno
Ncount = 1*rescount + A*AN + Cy*CyN + D*DN + E*EN + F*FN + G*GN + H*HN + I*IN + K*KN + L*LN + M*MN + N*NN + P*PN + Q*QN + R*RN + S*SN + T*TN + V*VN + W*WN + Y*YN
Ocount = 1 + 1*rescount + A*AO + Cy*CyO + D*DO + E*EO + F*FO + G*GO + H*HO + I*IO + K*KO + L*LO + M*MO + N*NO + P*PO + Q*QO + R*RO + S*SO + T*TO + V*VO + W*WO + Y*YO
Scount = A*AS + Cy*CyS + D*DS + E*ES + F*FS + G*GS + H*HS + I*IS + K*KS + L*LS + M*MS + N*NS + P*PS + Q*QS + R*RS + S*SS + T*TS + V*VS + W*WS + Y*YS
Htocount = Hexcount + Hnocount

MWul = Htocount*mwH + Ccount*mwC + Ncount*mwN + Ocount*mwO + Scount*mwS
MW2h = Hexcount*mwH + Hnocount*2 + Ccount*mwC + Ncount*mwN + Ocount*mwO + Scount*mwS
MW13c = Htocount*mwH + Ccount*13 + Ncount*mwN + Ocount*mwO + Scount*mwS
MW15n = Htocount*mwH + Ccount*mwC + Ncount*15 + Ocount*mwO + Scount*mwS
MW2h13c = Hexcount*mwH + Hnocount*2 + Ccount*13 + Ncount*mwN + Ocount*mwO + Scount*mwS
MW2h15n = Hexcount*mwH + Hnocount*2 + Ccount*mwC + Ncount*15 + Ocount*mwO + Scount*mwS
MW13c15n = Htocount*mwH + Ccount*13 + Ncount*15 + Ocount*mwO + Scount*mwS
MW2h13c15n = Hexcount*mwH + Hnocount*2 + Ccount*13 + Ncount*15 + Ocount*mwO + Scount*mwS

Abs280 = W*W280 + Y*Y280 
Abs280CyCy = Abs280 + CyCy*CyCy280 

Abs205 = BB*BB205 + W*W205 + F*F205 + Y*Y205 + H*H205 + M*M205 + R*R205 + N*N205 + Q*Q205 + Cy*Cy205
Abs205CyCy = Abs205 + CyCy*CyCy205

###############################################################################
# Output
###############################################################################

print("Content-type: text/plain")
print()
print("This script calculates molar absorptivities (extinction coefficients) at 205")
print("nm and 280 nm from an amino acid sequence. It also calculates the molecular")
print("weight for various universal isotopic labeling schemes.")
print()
print("Reference:")
print("Anthis N.J., Clore G.M. (2013) Sequence-specific determination of protein")
print("and peptide concentrations by absorbance at 205 nm, Protein Science 22, 851-8, doi:10.1002/pro.2253")
print()
print("Amino acid sequence:")
print(seq)
print()
print("Number of residues:")
print(rescount)
print()
print("Molecular weight (in H2O):")
print("natural abundance  %.2f" % (MWul))
print("2H                 %.2f" % (MW2h))
print("13C                %.2f" % (MW13c))
print("15N                %.2f" % (MW15n))
print("2H,13C             %.2f" % (MW2h13c))
print("2H,15N             %.2f" % (MW2h15n))
print("13C,15N            %.2f" % (MW13c15n))
print("2H,13C,15N         %.2f" % (MW2h13c15n))
print()
print("Amino acid composition:")
print("A = %s" % (A))
print("C = %s" % (Cy))
print("D = %s" % (D))
print("E = %s" % (E))
print("F = %s" % (F))
print("G = %s" % (G))
print("H = %s" % (H))
print("I = %s" % (I))
print("K = %s" % (K))
print("L = %s" % (L))
print("M = %s" % (M))
print("N = %s" % (N))
print("P = %s" % (P))
print("Q = %s" % (Q))
print("R = %s" % (R))
print("S = %s" % (S))
print("T = %s" % (T))
print("V = %s" % (V))
print("W = %s" % (W))
print("Y = %s" % (Y))
print()
print("Atomic composition:")
print("Carbon       (C)     %i" % (Ccount))
print("Hydrogen     (H)     %i" % (Htocount))
print("   non-exchangeable     %i" % (Hnocount))
print("   exchangeable         %i" % (Hexcount))
print("Nitrogen     (N)     %i" % (Ncount))
print("Oxygen       (O)     %i" % (Ocount))
print("Sulfur       (S)     %i" % (Scount))
print()
print("Molar absorptivity (extinction coefficient) at 280 nm =")
if CyCy >= 1:
 print("%i M-1 cm-1 (if no disulfide bonds are present)" % (Abs280))
 print("%i M-1 cm-1 (if all cysteines are in disulfide bonds)" % (Abs280CyCy))
else:
 print("%i M-1 cm-1" % (Abs280))
print()
print("Molar absorptivity (extinction coefficient) at 205 nm =")
if CyCy >= 1:
 print("%i M-1 cm-1 (if no disulfide bonds are present)" % (Abs205))
 print("%i M-1 cm-1 (if all cysteines are in disulfide bonds)" % (Abs205CyCy))
else:
 print("%i M-1 cm-1" % (Abs205))
print()
