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
    instances = odb.rootAssembly.instances['DREAM-1'.upper()]
    elementSets = instances.elementSets
    steps = odb.steps.keys()
    
    d = ScrapeODB(None, N, steps[0], odb=odb, instances=instances, elementSets=elementSets, frame_step=frame_step, num_grains=num_grains)
    for idx, step in enumerate(steps[1:]):
        print(step)
        try:
            di = ScrapeODB(None, N, step, odb=odb, instances=instances, elementSets=elementSets, frame_step=frame_step, num_grains=num_grains)
            d += di
        except TypeError:
            print('Exited on step {} (step {}/{}) - will attempt to save up to this step'.format(step, idx+2, len(steps)))
        
    return d
    

class ScrapeODB(): # Lattice
    def __init__(self, fpath, N=12, step='Loading',
                 frame_step=1, num_grains=None,
                 odb=None, instances=None, elementSets=None):
        """
        Args:
            fpath (str): Path to an odb file
            N (int): Number of slip systems
        """
        self.fpath = fpath
        if fpath == None:
            self.odb = odb
            self.instances = instances
            self.elementSets = elementSets
        else:
            self.odb = openOdb(path=self.fpath)    
            self.instances = self.odb.rootAssembly.instances['DREAM-1'.upper()] 
            self.elementSets = self.instances.elementSets
        #numElements=len(myInstance.elements) # Find all elements in grain
            
        self.step = step
        self.frames= self.odb.steps[step].frames
        self.num_frames = len(range(0, len(self.frames), frame_step))
        self.N = N
        
        ### Need to calc number of grains
        self.num_grains = len([i for i in self.elementSets.keys() if i[:5] == 'GRAIN'])
        self.num_grains = self.num_grains if num_grains == None else num_grains
        print(self.num_grains)
        
        # Create empy arrays to store data
        d_shape = (self.num_grains, self.num_frames)
        print(d_shape)
        
        sc = scrape_frames(self.frames, frame_step, self.num_grains, self.elementSets, self.N, self.step)
        self. s, self.e, self.lat, self.dims, self.rot, self.v, self.t = sc       

               
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


def scrape_frames(frames, frame_step, num_grains, elementSets, N, step):
    
    num_frames = len(frames)
    d_shape = (num_grains, num_frames)
    
    s = np.zeros((6,) + d_shape)
    e = np.zeros((6,) + d_shape)
    lat = np.zeros((6,) + d_shape)
    rot = np.zeros((3,) + d_shape)
    dims = np.zeros((3,) + d_shape)
    v = np.zeros(d_shape)
    t = np.zeros((num_frames,))
    
    lat_SDV_nums = range((N * 12 + 4), (N * 12 + 4) + 6 )
    rot_SDV_nums = range((N * 34 + 3), (N * 34 + 3) + 3 )
    
    for fidx in range(0, num_frames, frame_step): 
    
        frame = frames[fidx]
        data_fo = frame.fieldOutputs
    
        # Abaqus
        stress_ = data_fo['S']
        strain_ = data_fo['LE']
        Volume_ = data_fo['IVOL']
        
        co, la, ro = True, True, True
        try:
            coords_ = data_fo['COORD']
        except:
            co = False
        # SDVs
        try:
            latSDV_ = [data_fo['SDV{}'.format(i)] for i in lat_SDV_nums]
        except:
            la = False
        
        try:
            rotSDV_ = [data_fo['SDV{}'.format(i)] for i in rot_SDV_nums]
        except:
            ro = False
        
        for idx, grain in enumerate(range(1, num_grains + 1)):
            grain = 'GRAIN-{}'.format(grain)
            myInstance=elementSets[grain] # Select grain
            numElements=len(myInstance.elements) # Find all elements in grain
            
            # Elastic+plastic response
            stress = stress_.getSubset(region=myInstance,position=INTEGRATION_POINT,elementType='C3D8').values
            strain = strain_.getSubset(region=myInstance,position=INTEGRATION_POINT,elementType='C3D8').values
            Volume = Volume_.getSubset(region=myInstance,position=INTEGRATION_POINT,elementType='C3D8').values
                
            # SDV:  Lattice results exx, eyy, ezz, exy, exz, eyz
            if la:
                latSDV = [i.getSubset(region=myInstance,position=INTEGRATION_POINT,elementType='C3D8').values for i in latSDV_]
            if ro:
                rotSDV = [i.getSubset(region=myInstance,position=INTEGRATION_POINT,elementType='C3D8').values for i in rotSDV_]
            if co:
                coords = coords_.getSubset(region=myInstance,position=NODAL).values
                
    
            for ip in range(numElements*8):
                vv_ = Volume[ip].data
                
                s[:,idx, fidx] += stress[ip].data[:6] * vv_
                e[:,idx, fidx] += strain[ip].data[:6] * vv_
                
                if co:
                    coords[:,idx, fidx] += coords[ip].data[:3] * vv_
                if la:
                    lat[:,idx, fidx] += np.array([i[ip].data * vv_ for i in latSDV])
                if ro:
                    rot[:,idx, fidx] += np.array([i[ip].data * vv_ for i in rotSDV])
   
                v[idx, fidx] += vv_
    
            # Set values 
            vv = v[idx, fidx]
            s[:,idx, fidx] /= vv
            e[:,idx, fidx] /= vv
            lat[:,idx, fidx] /= vv       
            dims[:,idx, fidx] /= vv
            rot[:,idx, fidx] /= vv
            t[fidx] = frame.frameValue
            
            if fidx % 4 == 0:
                f = open("progress_{}.txt".format(num_grains),"w") 
                f.write('{}: {} out of {} frames complete'.format(step, fidx, num_frames)) 
                f.close()
            
    return  (s, e, lat, dims, rot, v, t)



t0 = time.time()

args.fpath = '/newhome/mi19356/chris_odb/chris_odb.odb' if args.fpath == None else args.fpath
args.N = 12 if args.N == None else int(args.N)
args.spath = "cpex_{}.npz".format(time.strftime("%Y%m%d_%H%M%S")) if args.spath == None else args.spath

if args.step == None:
    data = multi_step(args.fpath, args.N, num_grains=args.num_grains,
                     frame_step=args.frame_step)
else:
    data = ScrapeODB(args.fpath, args.N, args.step, num_grains=args.num_grains,
                     frame_step=args.frame_step)
    
try:
    data.save_cpex(args.spath)
except IOError:
    print('Invalid spath specified, saving as cpex_{data_time}.npz')
    data.save_cpex("cpex_{}.npz".format(time.strftime("%Y%m%d_%H%M%S")))

t1 = time.time()
np.savetxt('time_{}.txt'.format(data.num_grains), np.array([t1-t0]), delimiter=',')


