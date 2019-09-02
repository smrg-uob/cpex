import numpy as np
from argparse import ArgumentParser

parser = ArgumentParser()

parser.add_argument("-f", "--fpath", dest="fpath",
                    help="Path to file for scraping")

parser.add_argument("-N", "--nslip", dest="N",
                    help="Number of slip systems (default 12)")

parser.add_argument("-s", "--spath", dest="spath",
                    help="Path to save scraped data (default cpex.npz)")

args = parser.parse_args()

from odbAccess import *
from abaqusConstants import *
import time   




class ScrapeODB(): # Lattice
    def __init__(self, fpath, N=12):
        """
        Args:
            fpath (str): Path to an odb file
            N (int): Number of slip systems
        """
        self.fpath = fpath
        self.odb = openOdb(path=self.fpath)    
        self.frames= self.odb.steps['Loading'].frames
        self.instances = self.odb.rootAssembly.instances['DREAM-1'.upper()]
        print(len(self.frames))
        self.num_frames = len(self.frames)
        self.N = N
        
        ### Need to calc number of grains
        self.num_grains = 8 # HOW?

        # Create empy arrays to store data
        d_shape = (self.num_grains, self.num_frames)
        
        self.s = np.zeros((6,) + d_shape)
        self.e = np.zeros((6,) + d_shape)
        self.lat = np.zeros((6,) + d_shape)
        self.rot = np.zeros((3,) + d_shape)
        self.dims = np.zeros((3,) + d_shape)
        self.v = np.zeros(d_shape)

        # Open each frame in turn and extract data
        for fidx in range(self.num_frames): #self.num_frames): 
            frame = self.frames[fidx]
            sc = scrape_frame(frame, self.num_grains, self.instances, self.N)

            self.s[:, :, fidx] = sc[0]
            self.e[:, :, fidx] = sc[1]
            self.lat[:, :, fidx] = sc[2]
            self.dims[:, :, fidx] = sc[3]
            self.rot[:, :, fidx] = sc[4]
            self.v[:, fidx] = sc[5]   
            
            if fidx % 4 == 0:
                f = open("progress.txt","w") 
                f.write('{} out of {} frames complete'.format(fidx, self.num_frames)) 
                f.close()
               
            
    def save_cpex(self, fpath):
        np.savez(fpath, s=self.s, e=self.e, lat=self.lat, dims=self.dims, 
                 rot=self.rot, v=self.v, N=self.N,
                 num_frames=self.num_frames, num_grains=self.num_grains)


def scrape_frame(frame, num_grains, instances, N):
    
    lat_SDV_nums = range((N * 12 + 4), (N * 12 + 4) + 6 )
    rot_SDV_nums = range((N * 34 + 3), (N * 34 + 3) + 3 )
    
    s = np.zeros((6,num_grains))
    e = np.zeros((6,num_grains))  
    lat = np.zeros((6,num_grains))
    dims = np.zeros((3,num_grains))
    rot = np.zeros((3,num_grains))
    v = np.zeros((num_grains)) 
    
    data_fo = frame.fieldOutputs

    # Abaqus
    stress_ = data_fo['S']
    strain_ = data_fo['LE']
    Volume_ = data_fo['IVOL']
    coords_ = data_fo['COORD']
    
    # SDVs
    latSDV_ = [data_fo['SDV{}'.format(i)] for i in lat_SDV_nums]
    rotSDV_ = [data_fo['SDV{}'.format(i)] for i in rot_SDV_nums]
    
    for idx, grain in enumerate(range(1, num_grains + 1)):
        grain = 'GRAIN-{}'.format(grain)
        myInstance=instances.elementSets[grain] # Select grain
        numElements=len(myInstance.elements) # Find all elements in grain
        
        # Elastic+plastic response
        stress = stress_.getSubset(region=myInstance,position=INTEGRATION_POINT,elementType='C3D8').values
        strain = strain_.getSubset(region=myInstance,position=INTEGRATION_POINT,elementType='C3D8').values
        Volume = Volume_.getSubset(region=myInstance,position=INTEGRATION_POINT,elementType='C3D8').values
        coords = coords_.getSubset(region=myInstance,position=NODAL).values
        
        # SDV:  Lattice results exx, eyy, ezz, exy, exz, eyz
        latSDV = [i.getSubset(region=myInstance,position=INTEGRATION_POINT,elementType='C3D8').values for i in latSDV_]
        rotSDV = [i.getSubset(region=myInstance,position=INTEGRATION_POINT,elementType='C3D8').values for i in rotSDV_]
        
        # Iterate over total number of elements and sum the strain/unit volume
        s_v = np.zeros((6,))
        e_v = np.zeros((6,))
        lat_v = np.zeros((6,))
        coords_v = np.zeros((3,))
        rot_v = np.zeros((3,))
        vv = 0

        for ip in range(numElements*8):
            vv_ = Volume[ip].data
            
            s_v += np.array(stress[ip].data[:6] * vv_) 
            e_v += np.array(strain[ip].data[:6] * vv_) 
            
            # Handle this better!
            try:
                coords_v += np.array(coords[ip].data[:3] * vv_)
            except:
                pass
            
            try:
                lat_v += np.array([i[ip].data * vv_ for i in latSDV])
            except:
                pass
            
            try:
                rot_v += np.array([i[ip].data * vv_ for i in rotSDV])
            except:
                pass
            try:
                vv += vv_
            except:
                pass
        
        # Unpack values
        s[:,idx] = s_v / vv
        e[:,idx] = e_v  / vv
        lat[:,idx] = lat_v / vv            
        dims[:,idx] = coords_v / vv
        rot[:,idx] = rot_v / vv
        v[idx] = vv
        
    return  (s, e, lat, dims, rot, v)


t0 = time.time()

args.fpath = '/newhome/mi19356/chris_odb/chris_odb.odb' if args.fpath == None else args.fpath
args.N = 12 if args.N == None else args.N
args.spath = 'cpex.npz' if args.spath == None else args.spath

data = ScrapeODB(args.fpath, args.N)
data.save_cpex(args.spath)

t1 = time.time()
np.savetxt('time.txt', np.array([t1-t0]), delimiter=',')

