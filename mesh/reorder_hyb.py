###################################################################
# Copyright 2005 - 2021 Barcelona Supercomputing Center.          #
# Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE #
# for nonprofit scientific purposes only.                         #
# See companion file LICENSE.txt.                                 #
###################################################################



# Use after gmsh2alya for reordering HEX27 elements

import array as ar

oldfile="hex64.geo.dat"
newfile="replace.geo.dat"

f=open(oldfile,'r')
g=open(newfile,'w')

line = True
inELEM = False

while line:
      line = f.readline()
      if 'END_SKEW' in line:
         g.write(line)
         print('End Writting')
         print(line)
         line = False
      else:
         if 'ELEMENTS' in line:
            if 'END_ELEMENTS' in line:
               g.write(line)
               print('Finished element data!')
               inELEM = False
               line = True
            else:
               inELEM = True
               g.write(line)
               print('Reading element data...')
               line = True
         else:
            if inELEM:
               stripline = line.strip().split()
               print(stripline)
               s = [int(i) for i in stripline]
               print(s)
               nnodes = len(s)-1
               #print('nnodes = ', nnodes)

               if nnodes == 27:
                  #print('Reordering HEX27!')

                  a = [0] * 16
                  a[0] = s[12]
                  a[1] = s[14]
                  a[2] = s[10]
                  a[3] = s[17]
                  a[4] = s[19]
                  a[5] = s[20]
                  a[6] = s[18]
                  a[7] = s[11]
                  a[8] = s[13]
                  a[9] = s[15]
                  a[10] = s[16]
                  a[11] = s[23]
                  a[12] = s[24]
                  a[13] = s[22]
                  a[14] = s[25]
                  a[15] = s[21]

                  s[10] = a[0]
                  s[11] = a[1]
                  s[12] = a[2]
                  s[13] = a[3]
                  s[14] = a[4]
                  s[15] = a[5]
                  s[16] = a[6]
                  s[17] = a[7]
                  s[18] = a[8]
                  s[19] = a[9]
                  s[20] = a[10]
                  s[21] = a[11]
                  s[22] = a[12]
                  s[23] = a[13]
                  s[24] = a[14]
                  s[25] = a[15]

               elif nnodes == 18:
                  #print('Reordering PEN18!')

                  a = [0] * 7
                  a[0] = s[10]
                  a[1] = s[8]
                  a[2] = s[9]
                  a[3] = s[15]
                  a[4] = s[14]
                  a[5] = s[18]
                  a[6] = s[17]

                  s[8] = a[0]
                  s[9] = a[1]
                  s[10] = a[2]
                  s[14] = a[3]
                  s[15] = a[4]
                  s[17] = a[5]
                  s[18] = a[6]

               elif nnodes == 10:
                  #print('Reordering TET10!')

                  a = s[9]
                  b = s[10]
                  s[9] = b
                  s[10] = a

               elif nnodes == 14:
                  #print('Reordering PYR14!')

                  a = [0] * 5
                  a[0] = s[9]
                  a[1] = s[11]
                  a[2] = s[7]
                  a[3] = s[8]
                  a[4] = s[10]

                  s[7] = a[0]
                  s[8] = a[1]
                  s[9] = a[2]
                  s[10] = a[3]
                  s[11] = a[4]

               elif nnodes == 64:
                  print('Reordering HEX64!')

                  a = [0] * 48
                  a[0]  = s[10+1]
                  a[1]  = s[11+1]
                  a[2]  = s[13+1]
                  a[3]  = s[12+1]
                  a[4]  = s[15+1]
                  a[5]  = s[14+1]
                  a[6]  = s[16+1]
                  a[7]  = s[17+1]
                  a[8]  = s[18+1]
                  a[9]  = s[19+1]
                  a[10] = s[21+1]
                  a[11] = s[20+1]
                  a[12] = s[23+1]
                  a[13] = s[22+1]
                  a[14] = s[24+1]
                  a[15] = s[25+1]
                  a[16] = s[26+1]
                  a[17] = s[27+1]
                  a[18] = s[28+1]
                  a[19] = s[29+1]
                  a[20] = s[30+1]
                  a[21] = s[31+1]
                  a[22] = s[32+1]
                  a[23] = s[33+1]
                  a[24] = s[35+1]
                  a[25] = s[34+1]
                  a[26] = s[36+1]
                  a[27] = s[37+1]
                  a[28] = s[39+1]
                  a[29] = s[38+1]
                  a[30] = s[40+1]
                  a[31] = s[41+1]
                  a[32] = s[43+1]
                  a[33] = s[42+1]
                  a[34] = s[44+1]
                  a[35] = s[45+1]
                  a[36] = s[47+1]
                  a[37] = s[46+1]
                  a[38] = s[48+1]
                  a[39] = s[49+1]
                  a[40] = s[51+1]
                  a[41] = s[50+1]
                  a[42] = s[55+1]
                  a[43] = s[54+1]
                  a[44] = s[59+1]
                  a[45] = s[58+1]
                  a[46] = s[63+1]
                  a[47] = s[62+1]

                  s[14+1] = a[0]
                  s[15+1] = a[1]
                  s[18+1] = a[2]
                  s[19+1] = a[3]
                  s[11+1] = a[4]
                  s[10+1] = a[5]
                  s[24+1] = a[6]
                  s[25+1] = a[7]
                  s[28+1] = a[8]
                  s[29+1] = a[9]
                  s[30+1] = a[10]
                  s[31+1] = a[11]
                  s[27+1] = a[12]
                  s[26+1] = a[13]
                  s[12+1] = a[14]
                  s[13+1] = a[15]
                  s[16+1] = a[16]
                  s[17+1] = a[17]
                  s[20+1] = a[18]
                  s[21+1] = a[19]
                  s[22+1] = a[20]
                  s[23+1] = a[21]
                  s[40+1] = a[22]
                  s[43+1] = a[23]
                  s[42+1] = a[24]
                  s[41+1] = a[25]
                  s[44+1] = a[26]
                  s[45+1] = a[27]
                  s[46+1] = a[28]
                  s[47+1] = a[29]
                  s[36+1] = a[30]
                  s[37+1] = a[31]
                  s[38+1] = a[32]
                  s[39+1] = a[33]
                  s[49+1] = a[34]
                  s[48+1] = a[35]
                  s[51+1] = a[36]
                  s[50+1] = a[37]
                  s[32+1] = a[38]
                  s[35+1] = a[39]
                  s[34+1] = a[40]
                  s[33+1] = a[41]
                  s[54+1] = a[42]
                  s[55+1] = a[43]
                  s[58+1] = a[44]
                  s[59+1] = a[45]
                  s[62+1] = a[46]
                  s[63+1] = a[47]

               stripline = [str(i) for i in s]
               line = ' '.join(stripline)
               g.write(line+'\n')
               line = True
            else:
               g.write(line)
               line = True

f.close()
g.close()
