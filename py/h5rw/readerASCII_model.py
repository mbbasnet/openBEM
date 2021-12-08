'''
Created on Aug 4, 2017

@author: mbb
'''
import h5py
import numpy as np
import os
import sys

class ReaderASCII(object):
    '''
    classdocs
    '''

    def __init__(self, commentSign, commandSign):
        '''
        Constructor
        '''
        
        self.commentSign = commentSign
        self.commandSign = commandSign
        #
        self.currentLine = 0
        
        self.model_done = False # model command not done yet
        self.model_id = 0 # initiation of variable
        self.model_name = 'model'
        self.model_dim = 1
        
        self.nodeArray_label = []
        self.nodeArray_coord = []
        
        self.el_label = [] # element label
        self.element = []
        self.connectivity = []
        self.connectivityIndex = 1
        
        self.nodeset = []
        self.nodeset_name = []
        
        self.elset = []
        self.elset_name = []
        
        self.BEM_chunk = []
        self.material = []
        self.material_name = []

        self.step_type = []
        self.frequency = []
        
    def readInput(self, inputFile):
        '''
        reads the input file and prepares to be written in the HDF file
        '''
        self.inputFileName = inputFile
        # first filtering thie input (removing comments and blank lines, and saving as list)
        self.filterInput(inputFile) # returns refined input (in lowercase)
        ii = 0
        while ii < self.inputLength:
            rr = self.refined_input[ii]
            self.currentLine = ii
            print(ii, rr)
            if rr[0] == self.commandSign+'model':
                self.readModel()
            elif rr[0] == self.commandSign+'node':
                self.readNodes(self.model_dim)
                ii = self.currentLine # update the completed loop
            elif rr[0] == self.commandSign+'element':
                self.readElements()
                ii = self.currentLine # update the completed loop
            elif rr[0] == self.commandSign+'nset':
                self.readNodeset()
                ii = self.currentLine # update the completed loop
            elif rr[0] == self.commandSign+'elset':
                self.readElementset()
                ii = self.currentLine # update the completed loop
            elif rr[0] == self.commandSign+'bemchunk':
                self.readBEM_chunk()
                ii = self.currentLine # update the completed loop
            elif rr[0] == self.commandSign+'material':
                self.readMaterial()
                ii = self.currentLine # update the completed loop
            elif rr[0] == self.commandSign+'step':
                self.readStep()
                ii = self.currentLine # update the completed loop
            ii = ii+1
        print('Read input file compelted.')
                
    def filterInput(self, inputFile):
        '''
        inputFile = full path of input file
        comment lines and empty lines will be removed
        '''
        self.inputFile = inputFile # carry the file name and location
        inp = open(inputFile, 'r') # open the input file
        inpText = inp.read()  
        rows = inpText.split('\n') # separating the rows
        inp.close()
       
        reduced_rows = [] #inputFile[0:-3]]# first row will be hdf5 file name
        
        for rr in rows:
            rr=rr.replace(' ', '') # removing the whitespaces
            rr = rr.lower() # making all lowercase letters
            if rr != '':
                if rr[0] != self.commentSign:
                    rr = rr.split(',')
                    reduced_rows.append(rr)
        del rows
        self.refined_input = reduced_rows
        self.inputLength = len(reduced_rows)
        
    def readModel(self):
        '''
        Reads model entry from the input file
        '''
        rr = self.refined_input[self.currentLine]
        for kk in rr:
            kk = kk.split('=')
            if kk[0] == 'id': # the parameter is model ID
                try:
                    self.model_id = np.int32(kk[1])
                except ValueError:
                    print('Value given for ID is not an integer')
                    raise
            elif kk[0] == 'name':
                self.model_name = kk[1]
            elif kk[0] == 'dim':
                try:
                    modelDim = np.int32(kk[1])
                    self.model_dim =  modelDim 
                    # if model dimension not given, it should be read from number of coordinates in the nodes later
                except ValueError:
                    print('Value given for DIM is not an integer')
                    raise 
        print(self.model_id, self.model_name, self.model_dim)
        
    def readNodes(self, modelDim): #, line, label_array, coord_array):
        '''
        Reads nodes starting from the given line of the input file and saves to corresponding Arrays
        line: the index of the row in processed input file 
        label_array: array to store labels of the node number
        coord_array: array to store coordinates of the corresponding nodes
        '''
        
        # reading the node lines until the next command line
        jj = self.currentLine+1 # index of current line (next data line to the command node)
        nodeData = True # after starting of node command it is true that next data are node entries
        #
        while nodeData == True:
            dd = self.refined_input[jj] # the current line 
            ddTemp = dd[0] # to examine the first variable
            if ddTemp[0] == self.commandSign: # checking first letter of first entry in the row, new command starts
                print(ddTemp, 'in ', jj)
                nodeData = False # new command started
              
            else:
                try:
                    print(jj, 'node: ', self.refined_input[jj])
                    self.nodeArray_label.append(np.int32(dd[0]))
                    coord = []
                        
                    for kk in range(1,modelDim+1):
                        coord.append(np.float64(dd[kk]))
                    self.nodeArray_coord.append(coord)
                            
                except ValueError:
                    print('Value given for DIM is not a Real Value')
                    raise
                jj = jj+1
                if jj > len(self.refined_input)-1:
                    nodeData = False
            
        self.currentLine = jj-1 
        
    def readElements(self):
        '''
        reads each element in the row and saves to element array
        '''
        ee = Element()
        jj = self.currentLine # in the same line as element
        # now reading the attributes for the model
        rr = self.refined_input[jj]
        for kk in rr:
            kk = kk.split('=')
            if kk[0] == 'eltype' or kk[0] == 'type':
                try:
                    eltype_name = kk[1]
                    eltype_index, num_nodes = ee.element_ref(eltype_name) 
                    # gives corresponding integer for element type reference and number of nodes
                
                except ValueError:
                    print('the name is not element type name')
                    raise

        elementData = True
        jj = jj+1 # now start with data lines
        #

        while elementData == True: 
            dd = self.refined_input[jj] # the current line 
            ddTemp = dd[0] # to examine the first variable
            if ddTemp[0] == self.commandSign: # checking first letter of first entry in the row, new command starts
                elementData = False # new command started
            else:
                try:
                    
                    print(jj, 'element ')#,  self.refined_input[jj])
 
                    # [element label, eltype_reference, connectivity_start_index, connectivity_end_index]
                    self.el_label.append(np.int32(dd[0]))
                    self.element.append([np.int32(dd[0]), np.int32(eltype_index), \
                                         np.int32(self.connectivityIndex), np.int32(self.connectivityIndex+num_nodes-1)])
        
                
                    self.connectivityIndex = self.connectivityIndex+num_nodes
                    connectivity= []
                        
                    for kk in range(1,num_nodes+1):
                        # changing node labels to corresponding node index
                        node_label = np.int32(dd[kk])
                        #print('node label ', node_label, self.nodeArray_label)
                        node_index = self.nodeArray_label.index(node_label) + 1 # index based on 1 for fortran use
                        #print('node index ', node_index)
                        connectivity.append(np.int32(node_index))
                    # the collection of nodes as one dimensional arrays later to be reshaped in fortran
                    self.connectivity.extend(connectivity) # it is an one dimensional array
                    jj = jj+1
                except ValueError:
                    print('the name is not element type name')
                    raise
        self.currentLine = jj-1
        
        
    def readNodeset(self): #, line, label_array, coord_array):
        '''
        Reads nodes starting from the given line of the input file and saves to corresponding Arrays
        line: the index of the row in processed input file 
        label_array: array to store labels of the node number
        coord_array: array to store coordinates of the corresponding nodes
        '''
        # reading name of the nodeset
        rr = self.refined_input[self.currentLine]
        for ii in range(1, len(rr)):
            kk = rr[ii].split('=')
            if kk[0] == 'name' or kk[0] == 'nset':
                nset_name = kk[1]

        generate = False
        for ii in range(1, len(rr)):
            if rr[ii] == 'generate':
                generate = True
            

        # reading the node lines until the next command line
        jj = self.currentLine+1 # index of current line (next data line to the command node)
        nodesetData = True # after starting of node command it is true that next data are node entries
        nodeset = []

        if generate == True:
            try:
                ndata= self.refined_input[jj] # data line for nodeset
                nset_data = np.arange(np.int32(ndata[0]), np.int32(ndata[1])+1, np.int32(ndata[2])) # data of labels
                # now to get data of node indices
               
                for nn in nset_data:
                    node_index = self.nodeArray_label.index(np.int32(nn))+1
                    nodeset.append(np.int32(node_index))
                
            except ValueError:
                print('Error: density parameter not realistic')
                raise
            jj = jj+1

        elif generate == False:
            #
            while nodesetData == True:
                dd = self.refined_input[jj] # the current line 
                ddTemp = dd[0] # to examine the first variable
                if ddTemp[0] == self.commandSign: # checking first letter of first entry in the row, new command starts
                    print(ddTemp, 'in ', jj)
                    nodesetData = False # new command started
              
                else:
                    try:
                        #print(jj, 'nodeset: ', self.refined_input[jj])
                        for kk in dd:
                            node_index = self.nodeArray_label.index(np.int32(kk)) + 1 # node label changed to index starting from 1
                            nodeset.append(np.int32(node_index))
                    
                            
                    except ValueError:
                        print('Error while reading nodesets')
                        raise
                    jj = jj+1
                    if jj > len(self.refined_input)-1:
                        nodesetData = False

        self.nodeset.append(nodeset) # adding nodesets to the global list
        self.nodeset_name.append(nset_name) # appending the name of the node set
        print('the nodesets are:')
        print(self.nodeset_name)
        print(self.nodeset)
        self.currentLine = jj-1 
        
    def readElementset(self): #, line, label_array, coord_array):
        '''
        Reads nodes starting from the given line of the input file and saves to corresponding Arrays
        line: the index of the row in processed input file 
        label_array: array to store labels of the node number
        coord_array: array to store coordinates of the corresponding nodes
        '''
        # reading name of the nodeset
        rr = self.refined_input[self.currentLine]
        for ii in range(1, len(rr)):
            kk = rr[ii].split('=')
            if kk[0] == 'name' or kk[0] == 'elset':
                elset_name = kk[1]
        # reading the node lines until the next command line
        jj = self.currentLine+1 # index of current line (next data line to the command node)
        elsetData = True # after starting of node command it is true that next data are node entries
        elset = []
        #
        generate = False
        for ii in range(1, len(rr)):
            if rr[ii] == 'generate':
                generate = True
            

        if generate == True:
            try:
                eldata= self.refined_input[jj] # data line for nodeset
                elset_data = np.arange(np.int32(eldata[0]), np.int32(eldata[1])+1, np.int32(eldata[2])) # data of labels
                # now to get data of node indices
                elset = []
                for nn in elset_data:
                    el_index = self.el_label.index(np.int32(nn))+1
                    elset.append(np.int32(el_index))

            except ValueError:
                print('Error: density parameter not realistic')
                raise
            jj=jj+1
            
        elif generate == False:
            #
            while elsetData == True:
                dd = self.refined_input[jj] # the current line 
                ddTemp = dd[0] # to examine the first variable
                if ddTemp[0] == self.commandSign: # checking first letter of first entry in the row, new command starts
                    print(ddTemp, 'in ', jj)
                    elsetData = False # new command started
              
                else:
                    try:
                        print(jj, 'elset: ', self.refined_input[jj])
                        for kk in dd:
                            el_index = self.el_label.index(np.int32(kk)) + 1 # node label changed to index starting from 1
                            elset.append(np.int32(el_index))
                    
                            
                    except ValueError:
                        print('Error while reading nodesets')
                        raise
                    jj = jj+1
                    if jj > len(self.refined_input)-1:
                        elsetData = False

        self.elset.append(elset) # adding element sets to the global list
        self.elset_name.append(elset_name) # appending the name of the element set
        print('the elementsets are:')
        print(self.elset_name)
        print(self.elset)
        self.currentLine = jj-1 
        
    def readMaterial(self): #, line, label_array, coord_array):
        '''
        Reads nodes starting from the given line of the input file and saves to corresponding Arrays
        line: the index of the row in processed input file 
        label_array: array to store labels of the node number
        coord_array: array to store coordinates of the corresponding nodes
        '''
        # reading name of the nodeset
        mtrl = Material()
        MaterialData = True
        rr = self.refined_input[self.currentLine]
        for ii in range(1, len(rr)):
            kk = rr[ii].split('=')
            if kk[0] == 'name' or kk[0] == 'material':
                mtrl.name = kk[1]
                
        jj = self.currentLine+1 # go to sub commands
        while MaterialData == True:
            # reading the sub commands for material input
            rr = self.refined_input[jj] # new line subcommands of material
            
            if rr[0] == self.commandSign+'density':
                try:
                    jj = jj+1
                    data = self.refined_input[jj] # data line for density
                    mtrl.density = np.float64(data[0]) # only one input (density)
                    
                except ValueError:
                    print('Error: density parameter not realistic')
                    raise
                jj=jj+1
                
            elif rr[0] == self.commandSign+'elastic':
                try:
                    jj = jj+1
                    data = self.refined_input[jj] # data line for elastic material
                    mtrl.E = np.float64(data[0]) # Young's modulus of elasticity
                    mtrl.nu = np.float64(data[1]) # poisson's ratio
                    
                except ValueError:
                    print('Error: elastic parameters not realistic')
                    raise
                jj=jj+1
                
            elif rr[0][0]== self.commandSign:
                if rr[0][1:] not in mtrl.inputParam:
                    MaterialData = False
            else:
                MaterialData = False
            
        self.material.append(mtrl) # adding material to the material list
        self.material_name.append(mtrl.name)
        self.currentLine = jj-1 
        
                    
    def readBEM_chunk(self):
        '''
        Reads model entry from the input file
        '''
        bb = BEM_chunk()
        rr = self.refined_input[self.currentLine]
        for kk in rr:
            kk = kk.split('=')
            if kk[0] == 'bound': # the parameter is model ID
                try:
                    bb.bound = bb.domainType_ref(kk[1])
                except ValueError:
                    print('Value given is not listed in domain type dictionary')
                    raise
            elif kk[0] == 'element':
                # getting index of node set
                try:
                    elset_index = self.elset_name.index(kk[1])
                    bb.belm = self.elset[elset_index]
                except ValueError:
                    print('The element set is not found')
                    raise
            elif kk[0] == 'bnode':
                # getting index of node set
                try:
                    nset_index = self.nodeset_name.index(kk[1])
                    bb.bnode = self.nodeset[nset_index]
                except ValueError:
                    print('The nodeset is not found')
                    raise
            elif kk[0] == 'inode':
                # getting index of node set
                try:
                    nset_index = self.nodeset_name.index(kk[1])
                    bb.inode = self.nodeset[nset_index]
                except ValueError:
                    print('The nodeset is not found')
                    raise
            elif kk[0] == 'material':
                # getting index of node set
                try:
                    mat_index = self.material_name.index(kk[1])
                    bb.material = np.int32(mat_index+1) # index of the material starting from 1
                except ValueError:
                    print('The material does not exist')
                    raise
                
        self.BEM_chunk.append(bb)

    def readStep(self):
        '''
        Reads model entry from the input file
        '''
        jj = self.currentLine
        rr = self.refined_input[jj]
        for kk in rr:
            kk = kk.split('=')
            if kk[0] == 'type': # the parameter is model ID
                try:
                    if (kk[1] == 'frequency1' or (kk[1] == 'frequency3' or kk[1] == 'frequency2')):
                        jj = jj+1
                        rr = self.refined_input[jj]
                        if (len(rr) == 1 or np.int32(rr[2]) == 1):
                            self.frequency = [np.float32(rr[1])]
                        else:
                            intrvl = (np.float32(rr[1]) - np.float32(rr[0]))/(np.int32(rr[2])-1)
                            self.frequency = np.arange(np.float32(rr[0]), np.float32(rr[1])+intrvl/2.0, intrvl)
                        
                        for ll in range(0, len(self.frequency)):
                            if kk[1] == 'frequency1':
                                self.step_type.append(1) # 1 for antiplane case
                            elif kk[1] == 'frequency2':
                                self.step_type.append(2) # 2 for 2D plane case
                            elif kk[1] == 'frequency3':
                                self.step_type.append(3) # 3 for 3D case

                except ValueError:
                    print('Value given is not listed in domain type dictionary')
                    raise
        self.currentline = jj-1

                    
class CreaterHDF5(object):  
    '''
    creates HDF5 inputs
    '''
         
    def __init__(self):
        '''
        Constructor
        '''
        # the model attributes initiation
        self.model_id = 0
        
        # initiation of node arrays
        self.nodeArray_coord = []
        self.nodeArray_label = []
        
    def createModel_h5(self, inputModel, hdfFileName):
        '''
        input is instance of the the class readerASCII which will have the information required
        '''
        # first create HDF5 model
        fileName = './'+hdfFileName+'.h5'
        # check if file exists
        if os.path.isfile(fileName):
            file_overwrite = input('The file already exists. Overwrite? (Y/N)')
            if (file_overwrite.lower() == 'y' or file_overwrite.lower() =='yes'):
                os.remove(fileName)
            else:
                sys.exit('File already exists. New model not created')
               
        mdl = h5py.File(fileName,'a')
            
        # creating attributes in the model 
        mdl.attrs['id'] = np.int32(inputModel.model_id)
        mdl.attrs['name'] = inputModel.model_name
        mdl.attrs['dim'] = inputModel.model_dim
        mdl.attrs['num_material'] = np.int32(len(inputModel.material_name))
        mdl.attrs['num_BEM_chunk'] = np.int32(len(inputModel.BEM_chunk))
        mdl.attrs['num_step'] = np.int32(len(inputModel.step_type))
        
        # creating mesh
        mesh = mdl.create_group('mesh') # creating mesh group
        node = mesh.create_group('node') # creating node group
        elm = mesh.create_group('element') # creating node group
        nset = mesh.create_group('nodeset') 
        elset = mesh.create_group('elementset')

        
        # creating nodes and attributes
        mesh.attrs['dim'] = np.int32(len(inputModel.nodeArray_coord[0]))
        mesh.attrs['num_node'] = np.int32(len(inputModel.nodeArray_label))
        mesh.attrs['num_nodeset'] = np.int32(len(inputModel.nodeset))
        mesh.attrs['num_elementset'] = np.int32(len(inputModel.elset))
       
        node.create_dataset('label', shape=(np.int32(len(inputModel.nodeArray_label)),1), \
                            dtype= np.int32, data= inputModel.nodeArray_label)
        node.create_dataset('coord', shape=(np.int32(len(inputModel.nodeArray_coord)),np.int32(len(inputModel.nodeArray_coord[0]))), \
                            dtype= np.float64, data= inputModel.nodeArray_coord)
        
        # creating elements and attributes
        mesh.attrs['num_element'] =  np.int32(len(inputModel.element))
        mesh.attrs['num_connectivity'] = np.int32(len(inputModel.connectivity))

        elm.create_dataset(name='index', shape=(len(inputModel.element), len(inputModel.element[0])),\
                            dtype= np.int32, data= inputModel.element)
        
        elm.create_dataset(name='connectivity', shape=(len(inputModel.connectivity), 1),\
                            dtype= np.int32, data= inputModel.connectivity)
        
        # assigning the nodesets
        if inputModel.nodeset:
            for ii in range(0,len(inputModel.nodeset)):
                dset_name = 'nset'+np.str(ii+1)
                dset = nset.create_dataset(name=dset_name, shape=(len(inputModel.nodeset[ii]), 1),\
                            dtype= np.int32, data= inputModel.nodeset[ii])
                dset.attrs['name'] = inputModel.nodeset_name[ii]
                dset.attrs['len'] = len(inputModel.nodeset[ii])
                
        # assigning the element sets
        if inputModel.elset:
            for ii in range(0,len(inputModel.elset)):
                dset_name = 'elset'+np.str(ii+1)
                dset = elset.create_dataset(name=dset_name, shape=(len(inputModel.elset[ii]), 1),\
                            dtype= np.int32, data= inputModel.elset[ii])
    
                dset.attrs['name'] = inputModel.elset_name[ii]
                dset.attrs['len'] = len(inputModel.elset[ii])
        
        # assigning BEM_chunks        
        if inputModel.BEM_chunk:
            bem_chunk = mdl.create_group('bem_chunk')
            for ii in range(0,len(inputModel.BEM_chunk)):
                chunk_name = 'bem'+np.str(ii+1)
                bchunk = bem_chunk.create_group(chunk_name) # creating a subgroup
                # now writing data to bem chunks
                bchunk.attrs['bound'] = np.int32(inputModel.BEM_chunk[ii].bound)
                bchunk.attrs['material'] = np.int32(inputModel.BEM_chunk[ii].material)
                bchunk.attrs['num_bnode'] = np.int32(len(inputModel.BEM_chunk[ii].bnode))
                bchunk.attrs['num_inode'] = np.int32(len(inputModel.BEM_chunk[ii].inode))
                bchunk.attrs['num_belm'] = np.int32(len(inputModel.BEM_chunk[ii].belm))
                
                bchunk.create_dataset(name='bnode', shape=(len(inputModel.BEM_chunk[ii].bnode),1),\
                                       dtype=np.int32, data=inputModel.BEM_chunk[ii].bnode)
                
                if bchunk.attrs['num_inode'] != np.int32(0):
                    bchunk.create_dataset(name='inode', shape=(len(inputModel.BEM_chunk[ii].inode),1),\
                                       dtype=np.int32, data=inputModel.BEM_chunk[ii].inode)
                
                bchunk.create_dataset(name='belm', shape=(len(inputModel.BEM_chunk[ii].belm),1),\
                                       dtype=np.int32, data=inputModel.BEM_chunk[ii].belm)
                
        # assigning the materials
        if inputModel.material:
            mtrl_group = mdl.create_group('material')
            for ii in range(0,len(inputModel.material_name)):
                mat_name = 'mat'+np.str(ii+1)
                mtrl = mtrl_group.create_group(mat_name)
                mtrl.attrs['name'] = inputModel.material_name[ii]
                if inputModel.material[ii].density:
                    mtrl.create_dataset(name='density', dtype=np.float64, data=inputModel.material[ii].density)
                    
                if inputModel.material[ii].E:
                    mtrl.create_dataset(name='E', dtype=np.float64, data=inputModel.material[ii].E)
                    
                if inputModel.material[ii].nu:
                    mtrl.create_dataset(name='nu', dtype=np.float64, data=inputModel.material[ii].nu)

        # assigning the materials
        if inputModel.step_type:
            step_group = mdl.create_group('step')
            for ii in range(0,len(inputModel.step_type)):
               step_name = 'step'+np.str(ii+1)
               step = step_group.create_group(step_name)
               step.create_dataset(name='step_type', dtype=np.int32, data=inputModel.step_type[ii])
               step.create_dataset(name='frequency', dtype=np.float64, data=inputModel.frequency[ii])

        print('hdf input created')
        
        
class Element(object):
    '''
    Referenced to element related activities
    '''
    
    def __init__(self):
        '''
        element types are referenced as integer
        '''
        
    
    def element_ref(self,el_name):
        '''
        corresponding integer id for element name
        '''
        # using element dictionary
        eltypeIndex = {'ps2dl0':1, 'ap2dl0': 2, '3dtry4': 3}
        eltypeNodes = {'ps2dl0':2, 'ap2dl0': 2, '3dtry4': 4}
        
        # returning element type index and number of nodes
        index = eltypeIndex[el_name]
        numNodes = eltypeNodes[el_name]
        
        return index, numNodes
        
        
class BEM_chunk(object):
    '''
    a BEM chunk object
    '''
    def __init__(self):
        '''
        '''
       
        self.bnode = [] # boundary nodes
        self.belm = [] # boundary elements
        self.inode = [] # internal nodes
        self.bound = 0
        self.material = 0
        
    def domainType_ref(self,bound):
        # dictionary for domain type
        domainType = {'undefined':0, 'finite':0, 'semiinfinite':1, 'infinite':2, 'smooth':4}
        
        bound_id = domainType[bound]
        
        return bound_id
    
class Material(object):
    '''
    material class
    '''
    def __init__(self):
        self.name = 'material'
        self.density = np.float64(0.)
        self.E = np.float64(0.)
        self.nu = np.float64(0.)
        
        self.define_inputParameters() # define set of input parameters for input file
        
    def define_inputParameters(self):
        self.inputParam = {'density', 'elastic'}