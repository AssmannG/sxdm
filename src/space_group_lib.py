'''
created by S.Basu
Date: 20-Oct-2017
'''

import os, sys
'''
a handy dictionary to map space group number and name along with lattice types
for 65 chiral space groups (proteins only..)
aP = triclinic Primitive,
mP = monoclinic Primitive, mC = monoclinic C-centered,
oP = orthorhombic Primitive, oC = orthorhombic C-centered,
oF = orthorhombic Face-centered, oI = orthorhombic body-centered
tP = tetragonal Primitive, tI = tetragonal body-centered
trP = trigonal Primitive
Rh = Rhombohedral hexagonal face
hP = hexagonal Primitive
cP = cubic Primitive, cF = cubic Face-centered, cI = cubic body-centered
'''
space_group = {
'1':('P1','aP'),
'3':('P121','mP'),
'4':('P1211','mP'),
'5':('C121','mC'),
'16':('P222','oP'),
'17':('P2221','oP'),
'18':('P21212','oP'),
'19':('P212121','oP'),
'20':('C2221','oC'),
'21':('C222', 'oC'),
'22':('F222', 'oF'),
'23':('I222', 'oI'),
'24':('I212121', 'oI'),
'75':('P4','tP'),
'76':('P41','tP'),
'77':('P42','tP'),
'78':('P43','tP'),
'79':('I4', 'tI'),
'80':('I41','tI'),
'89':('P422','tP'),
'90':('P4212','tP'),
'91':('P4122','tP'),
'92':('P41212','tP'),
'93':('P4222','tP'),
'94':('P42212','tP'),
'95':('P4322','tP'),
'96':('P43212','tP'),
'97':('I422','tI'),
'98':('I4122','tI'),
'143':('P3','trP'),
'144':('P31','trP'),
'145':('P32','trP'),
'146':('R3','Rh'),
'149':('P312','trP'),
'150':('P321','trP'),
'151':('P3112','trP'),
'152':('P3121','trP'),
'153':('P3212','trP'),
'154':('P3221','trP'),
'155':('R32','Rh'),
'168': ('P6','hP'),
'169': ('P61','hP'),
'170': ('P65','hP'),
'171': ('P62','hP'),
'172': ('P64','hP'),
'173': ('P63','hP'),
'177': ('P622','hP'),
'178':('P6122','hP'),
'179': ('P6522','hP'),
'180':('P6222','hP'),
'181':('P6422','hP'),
'182':('P6322','hP'),
'195':('P23','cP'),
'196':('F23','cF'),
'197':('I23', 'cI'),
'198':('P213','cP'),
'199':('I213', 'cI'),
'207':('P432', 'cP'),
'208':('P4232','cP'),
'209':('F432','cF'),
'210':('F4132', 'cF'),
'211':('I432', 'cI'),
'212':('P4332', 'cP'),
'213':('P4132','cP'),
'214':('I4132','cI')
}
