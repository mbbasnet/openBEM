'''
Created on Aug 4, 2017

@author: mbb
'''
import sys
sys.path.insert(0, "./py")
from h5rw.readerASCII_model import ReaderASCII,CreaterHDF5

if __name__ == '__main__':
    pass

x = 'mbb.inp'
print(x[-4:])
rd= ReaderASCII(commentSign='!', commandSign='*')
hd = CreaterHDF5()

# first reading input file
rd.readInput('./inp/firstTry.inp') # reading ASCII and putting into rows
hd.createModel_h5(rd, 'inp/firstTry')

