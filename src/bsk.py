import numpy as np
import numpy.linalg as la
from tqdm import tqdm

def NumberConnections(Max,a,b): #Compute the ID, number of connections, and the connections
    pbar = tqdm(total = Max, desc="Computing Connections")
    ID=np.zeros(Max)
    Nc=np.zeros(Max)
    Connections=[]
    for n in range(Max):
        kx=b[(a==n)]
        ky=a[(b==n)]
        ID[n]=n
        Nc[n]=len(kx)+len(ky)
        Connections.append(np.concatenate([kx,ky],axis=0))
        pbar.update(1)
    pbar.close()
    return ID.astype(int),Nc.astype(int),np.array(Connections)

def DP(x1,y1,z1,x2,y2,z2):  #Compute the distance between two points
    return np.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)

def Eigen(A):  # Return the eigenvalues of a Matrix
    w,v = la.eig(A)
    w[(w<0)]=0.0001   # To fix eigenvalues -10**-26
    return w

def Volume(a,b,c): # Compute the volume using the values, a,b,c
    v = a*b*c
    if v<=0:
        v=0.0001 # To fix densisty inf
    return v

def Features(Max,ID,x,y,z,nc,con): #Return the values of average distance and a,b,c values
    #------------ Nodes
    ad = np.zeros(Max)
    a = np.zeros(Max)
    b = np.zeros(Max)
    c = np.zeros(Max)
    vol = np.zeros(Max)
    den = np.zeros(Max)+np.exp(60)
    IDnn = ID[nc!=0]
    pbar = tqdm(total=len(IDnn), desc="Computing Features")
    for n in IDnn:
        Ixx = 0
        Iyy = 0
        Izz = 0
        Ixy = 0
        Iyz = 0
        Ixz = 0
        dist_temp = 0
        for i in con[n]:
            dist_temp = dist_temp+DP(x[n],y[n],z[n],x[i],y[i],z[i])
            Ixx = Ixx+(y[i]-y[n])**2+(z[i]-z[n])**2
            Iyy = Iyy+(z[i]-z[n])**2+(x[i]-x[n])**2
            Izz = Izz+(x[i]-x[n])**2+(y[i]-y[n])**2
            Ixy = Ixy-(x[i]-x[n])*(y[i]-y[n])
            Iyz = Iyz-(y[i]-y[n])*(z[i]-z[n])
            Ixz = Ixz-(x[i]-x[n])*(z[i]-z[n])
        pbar.update(1)
        
        ad[n] = dist_temp/nc[n]
        A=np.array([[Ixx,Ixy,Ixz],[Ixy,Iyy,Iyz],[Ixz,Iyz,Izz]])
        eig = Eigen(A)
        eig = np.sort(eig)[::-1]# Major -> Minor
        a[n] = np.sqrt(eig[0]).real #.real
        b[n] = np.sqrt(eig[1]).real
        c[n] = np.sqrt(eig[2]).real
        vol[n] = Volume(a[n],b[n],c[n])
        den[n] = 1/vol[n]
    pbar.close()
    return ad,vol,den

def neigh_features(Max,ID,nc,ncn,ad,den,con): #Compute all properties for first neighbors (Gradient)
    nc_n = np.zeros(Max)
    ad_n = np.zeros(Max)
    a_n = np.zeros(Max)
    b_n = np.zeros(Max)
    c_n = np.zeros(Max)
    den_n = np.zeros(Max)
    IDnn = ID[(nc!=0)]
    kk=0
    pbar = tqdm(total=len(IDnn), desc="Computing Delta Features")
    for i in IDnn:
        nc_temp = 0.0
        ad_temp = 0.0
        den_temp = 0.0
        for j in con[i]:
            nc_temp = nc_temp+ncn[j]
            ad_temp = ad_temp+ad[j]
            den_temp = den_temp+den[j]
        nc_n[i] = nc_temp/(1.0*nc[i])
        ad_n[i] = ad_temp/(1.0*nc[i])
        den_n[i] = den_temp/(1.0*nc[i])
        kk = kk+1
        pbar.update(1)     
    pbar.close()
    nc_n = nc_n-ncn
    ad_n = ad_n-ad
    den_n = den_n-den
    return nc_n,ad_n,den_n