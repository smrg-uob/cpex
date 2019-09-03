import numpy as np
from argparse import ArgumentParser

parser = ArgumentParser()

parser.add_argument("-f", "--fpath", dest="fpath",
                    help="Path to file for scraping")

parser.add_argument("-N", "--nslip", dest="N",
                    help="Number of slip systems (default 12)", type=int)

parser.add_argument("-s", "--spath", dest="spath",
                    help="Path to save scraped data (default cpex.npz)")

parser.add_argument("--step", dest="step",
                    help="Analysis step to extract from (default all steps)")

parser.add_argument("--fstep", dest="frame_step",
                    help="Step between frames (default 1)", type=int)

parser.add_argument("--ngrains", dest="num_grains",
                    help="Number of grains to analyse (default all)", type=int)

args = parser.parse_args()
args.frame_step = 1 if args.frame_step == None else args.frame_step


from odbAccess import *
from abaqusConstants import *
import time   


def multi_step(fpath, N, frame_step=1, num_grains=None):
    
    odb = openOdb(path=fpath)
    steps = odb.steps.keys()
    
    d = ScrapeODB(None, N, steps[0], odb, frame_step=frame_step, num_grains=num_grains)
    for idx, step in enumerate(steps[1:]):
        print(step)
        try:
            di = ScrapeODB(None, N, step, odb, frame_step=frame_step, num_grains=num_grains)
            d += di
        except TypeError:
            print('Exited on step {} (step {}/{}) - will attempt to save up to this step'.format(step, idx+2, len(steps)))
        
    return d
    

class ScrapeODB(): # Lattice
    def __init__(self, fpath, N=12, step='Loading', odb=None, 
                 frame_step=1, num_grains=None):
        """
        Args:
            fpath (str): Path to an odb file
            N (int): Number of slip systems
        """
        self.fpath = fpath
        if fpath == None:
            self.odb = odb
        else:
            self.odb = openOdb(path=self.fpath)    
        self.step = step
        
        self.frames= self.odb.steps[step].frames
        self.instances = self.odb.rootAssembly.instances['DREAM-1'.upper()]
        self.num_frames = len(self.frames) // frame_step
        self.N = N
        
        ### Need to calc number of grains
        self.num_grains = len([i for i in self.instances.elementSets.keys() if i[:5] == 'GRAIN'])
        self.num_grains = self.num_grains if num_grains == None else num_grains
        print(self.num_grains)
        
        # Create empy arrays to store data
        d_shape = (self.num_grains, self.num_frames)
        print(d_shape)
        
        self.s = np.zeros((6,) + d_shape)
        self.e = np.zeros((6,) + d_shape)
        self.lat = np.zeros((6,) + d_shape)
        self.rot = np.zeros((3,) + d_shape)
        self.dims = np.zeros((3,) + d_shape)
        self.v = np.zeros(d_shape)
        self.t = np.zeros((self.num_frames,))

        # Open each frame in turn and extract data
        for fidx in range(0, self.num_frames, frame_step): #self.num_frames): 
            frame = self.frames[fidx]
            sc = scrape_frame(frame, self.num_grains, self.instances, self.N)

            self.s[:, :, fidx] = sc[0]
            self.e[:, :, fidx] = sc[1]
            self.lat[:, :, fidx] = sc[2]
            self.dims[:, :, fidx] = sc[3]
            self.rot[:, :, fidx] = sc[4]
            self.v[:, fidx] = sc[5]  
            self.t[fidx] = sc[6]
            
            if fidx % 4 == 0:
                f = open("progress_10g.txt","w") 
                f.write('{}: {} out of {} frames complete'.format(self.step, fidx, self.num_frames)) 
                f.close()
        

               
            
    def save_cpex(self, fpath):
        np.savez(fpath, s=self.s, e=self.e, lat=self.lat, dims=self.dims, 
                 rot=self.rot, v=self.v, N=self.N,
                 num_frames=self.num_frames, time=self.t,
                 num_grains=self.num_grains)
        
    def __add__(self, other):
        
        self.num_frames += other.num_frames
        self.s = np.append(self.s, other.s, axis=-1)
        self.e = np.append(self.e, other.e, axis=-1)
        self.lat = np.append(self.lat, other.lat, axis=-1)
        self.dims = np.append(self.dims, other.dims, axis=-1)
        self.rot = np.append(self.rot, other.rot, axis=-1)
        self.v = np.append(self.v, other.v, axis=-1)
        self.t = np.append(self.t, other.t + self.t[-1] + other.t[1]/1e6, axis=-1)
        
        return self
        #return basic_merge([self, other])
        


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
    try:
        coords_ = data_fo['COORD']
    except:
        pass
    # SDVs
    try:
        latSDV_ = [data_fo['SDV{}'.format(i)] for i in lat_SDV_nums]
    except:
        pass
    
    try:
        rotSDV_ = [data_fo['SDV{}'.format(i)] for i in rot_SDV_nums]
    except:
        pass
    
    for idx, grain in enumerate(range(1, num_grains + 1)):
        grain = 'GRAIN-{}'.format(grain)
        myInstance=instances.elementSets[grain] # Select grain
        numElements=len(myInstance.elements) # Find all elements in grain
        
        # Elastic+plastic response
        stress = stress_.getSubset(region=myInstance,position=INTEGRATION_POINT,elementType='C3D8').values
        strain = strain_.getSubset(region=myInstance,position=INTEGRATION_POINT,elementType='C3D8').values
        Volume = Volume_.getSubset(region=myInstance,position=INTEGRATION_POINT,elementType='C3D8').values
        
        
        # SDV:  Lattice results exx, eyy, ezz, exy, exz, eyz
        try:
            latSDV = [i.getSubset(region=myInstance,position=INTEGRATION_POINT,elementType='C3D8').values for i in latSDV_]
        except:
            pass
        
        try:
            rotSDV = [i.getSubset(region=myInstance,position=INTEGRATION_POINT,elementType='C3D8').values for i in rotSDV_]
        except:
            pass
        
        try:
            coords = coords_.getSubset(region=myInstance,position=NODAL).values
        except:
            pass
            
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


            vv += vv_

        
        # Unpack values
        s[:,idx] = s_v / vv
        e[:,idx] = e_v  / vv
        lat[:,idx] = lat_v / vv            
        dims[:,idx] = coords_v / vv
        rot[:,idx] = rot_v / vv
        v[idx] = vv
        t = frame.frameValue
        
    return  (s, e, lat, dims, rot, v, t)



t0 = time.time()

args.fpath = '/newhome/mi19356/chris_odb/chris_odb.odb' if args.fpath == None else args.fpath
args.N = 12 if args.N == None else int(args.N)
args.spath = 'cpex.npz' if args.spath == None else args.spath

if args.step == None:
    data = multi_step(args.fpath, args.N, num_grains=args.num_grains,
                     frame_step=args.frame_step)
else:
    data = ScrapeODB(args.fpath, args.N, args.step, num_grains=args.num_grains,
                     frame_step=args.frame_step)
    
data.save_cpex(args.spath)

t1 = time.time()
np.savetxt('time_10g.txt', np.array([t1-t0]), delimiter=',')

