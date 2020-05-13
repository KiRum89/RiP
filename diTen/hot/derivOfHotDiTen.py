import numpy as np
from hotDiTen import HotDiTen
from selector import H
from scipy.special import iv, ivp



class DerivOfHotDiTen_num:

   

    def __init__(self,X,Yabs,gamma,local_Nz,local_Nx,num_term):
        
        self.X = X
        self.Yabs = Yabs
        self.gamma = gamma
        self.Nz = local_Nz
        self.Nx = local_Nx
        self.num_term = num_term

        
        
    def d_dX(self):
        hotDiTen = ["Kxx","Cxy","Cxz","Kyy","Cyz","Kzz"]
        dX = 1e-6
        val = np.array([])
        for el in hotDiTen:
            val1=HotDiTen.__dict__[el]( HotDiTen(self.X,self.Yabs,self.gamma,self.Nz,self.Nx,self.num_term) )
            val2=HotDiTen.__dict__[el]( HotDiTen(self.X+dX,self.Yabs,self.gamma,self.Nz,self.Nx,self.num_term) )
            
            val = np.concatenate((val,val2-val1))
    
        return val/dX

    def d_dY(self):
        hotDiTen = ["Kxx","Cxy","Cxz","Kyy","Cyz","Kzz"]
        dY = 1e-6
        val = np.array([])
        for el in hotDiTen:
            val1=HotDiTen.__dict__[el]( HotDiTen(self.X,self.Yabs,self.gamma,self.Nz,self.Nx,self.num_term) )
            val2=HotDiTen.__dict__[el]( HotDiTen(self.X,self.Yabs+dY,self.gamma,self.Nz,self.Nx,self.num_term) )
            

            val = np.concatenate((val,val2-val1))
    
        return val/dY


    def d_dgamma(self):
        hotDiTen = ["Kxx","Cxy","Cxz","Kyy","Cyz","Kzz"]
        dgamma = 1e-6
        val = np.array([])
        for el in hotDiTen:
            val1=HotDiTen.__dict__[el]( HotDiTen(self.X,self.Yabs,self.gamma,self.Nz,self.Nx,self.num_term) )
            val2=HotDiTen.__dict__[el]( HotDiTen(self.X,self.Yabs,self.gamma+dgamma,self.Nz,self.Nx,self.num_term) )
            
            val = np.concatenate((val,val2-val1))
    
        return val/dgamma


    def d_dNz(self):
        hotDiTen = ["Kxx","Cxy","Cxz","Kyy","Cyz","Kzz"]
        dNz = 1e-3
        val = np.array([])
        for el in hotDiTen:
            val1=HotDiTen.__dict__[el]( HotDiTen(self.X,self.Yabs,self.gamma,self.Nz,self.Nx,self.num_term) )
            val2=HotDiTen.__dict__[el]( HotDiTen(self.X,self.Yabs,self.gamma,self.Nz+dNz,self.Nx,self.num_term) )
            
            val = np.concatenate((val,val2-val1))
    
        return val/dNz



    def d_dNx(self):
        hotDiTen = ["Kxx","Cxy","Cxz","Kyy","Cyz","Kzz"]
        dNx = 1e-6
        val = np.array([])
        for el in hotDiTen:
            val1=HotDiTen.__dict__[el]( HotDiTen(self.X,self.Yabs,self.gamma,self.Nz,self.Nx,self.num_term) )
            val2=HotDiTen.__dict__[el]( HotDiTen(self.X,self.Yabs,self.gamma,self.Nz,self.Nx+dNx,self.num_term) )
            
            val = np.concatenate((val,val2-val1))
    
        return val/dNx


## analitical derivatives of diel tensor            
class DerivOfHotDiTen(HotDiTen):


    def __init__(self,X,Yabs,gamma,local_Nz,local_Nx,num_term):

        self.X = X
        self.Yabs = Yabs
        self.gamma = gamma
        self.Nz = local_Nz
        self.Nx = local_Nx
        self.num_term = num_term
        self.zeta0 = self.zeta(0)
        self.n= np.arange(1,num_term+1)
        self.l = self.lam()
        


    def dlam_du(self,u):
        return 2*self.l*(1/self.Nx*H(u,"Nx") + 1/self.gamma*H(u,"gamma") - 1/self.Yabs*H(u,"Y"))

        
    def dzetan_du(self,u,n):
        zeta_n = self.zeta(n) 
        #print (1+n*self.Yabs)
        #if (1+n*self.Yabs)
        return zeta_n*(n/(1+n*self.Yabs)*H(u,"Y")  - 1/self.Nz*H(u,"Nz") - 1/self.gamma*H(u,"gamma")  )

    def dKxx_du(self,u):

        sum=np.sum(self.n**2*( self.dlam_du(u)*ivp(self.n,self.l)*(self.Z(self.n) + self.Z(-self.n)) + iv(self.n,self.l)*(self.dzetan_du(u,self.n)*self.Zp(self.n)+self.dzetan_du(u,-self.n)*self.Zp(-self.n))))
        return self.pdBxx_du(u)*(self.Kxx()-1) + self.Bxx()*sum
    

    def dCxy_du(self,u):
        sum=np.sum(self.n*( self.dlam_du(u)*(ivp(self.n,self.l,1) - ivp(self.n,self.l,2))*(self.Z(self.n)-self.Z(-self.n))+(self.dzetan_du(u,self.n)*self.Zp(self.n)-self.dzetan_du(u,-self.n)*self.Zp(-self.n))*(iv(self.n,self.l)-ivp(self.n,self.l))))

        
        return self.pdBxy_du(u)*self.Cxy() - self.Bxy()*sum

    
    def dCxz_du(self,u):

    
        sum=np.sum(self.n*( self.dlam_du(u)*ivp(self.n,self.l,1)*(-self.Zp(self.n)+self.Zp(-self.n))+iv(self.n,self.l)*(-self.dzetan_du(u,self.n)*self.Zpp(self.n)+self.dzetan_du(u,-self.n)*self.Zpp(-self.n)))/2)

        return self.pdBxz_du(u)*self.Cxz() - self.Bxz()*sum




    def dKyy_du(self,u):

        
    
        sum= np.sum(self.dlam_du(u) * ( self.n**2*ivp(self.n,self.l) + 4*self.l*( iv(self.n,self.l)-ivp(self.n,self.l) ) + 2*self.l**2*(ivp(self.n,self.l)-ivp(self.n,self.l,2)) )*(self.Z(self.n)+self.Z(-self.n)) + (self.dzetan_du(u,self.n)*self.Zp(self.n)+self.dzetan_du(u,-self.n)*self.Zp(-self.n))*(self.n**2*iv(self.n,self.l)+2*self.l**2*(iv(self.n,self.l)-ivp(self.n,self.l))))
        
        nn = 0
        firstTerm=self.dlam_du(u) * ( nn**2*ivp(nn,self.l) + 4*self.l*( iv(nn,self.l)-ivp(nn,self.l) ) + 2*self.l**2*(ivp(nn,self.l)-ivp(nn,self.l,2) ))*(self.Z(nn)+self.Z(-nn)) + (self.dzetan_du(u,nn)*self.Zp(nn)+self.dzetan_du(u,-nn)*self.Zp(-nn))*(nn**2*iv(nn,self.l)+2*self.l**2*(iv(nn,self.l)-ivp(nn,self.l)))

        return self.pdByy_du(u)*(self.Kyy()-1) + self.Byy()*(sum + 0.5*firstTerm)


    def dCyz_du(self,u):
    

        sum=np.sum(self.dlam_du(u)*(ivp(self.n,self.l)-ivp(self.n,self.l,2))*(self.zeta(self.n)*self.Z(self.n)+self.zeta(-self.n)*self.Z(-self.n))\
               +(iv(self.n,self.l)-ivp(self.n,self.l))*(self.dzetan_du(u,self.n)*(self.Z(self.n)+self.zeta(self.n)*self.Zp(self.n)) + self.dzetan_du(u,-self.n)*(self.Z(-self.n)+self.zeta(-self.n)*self.Zp(-self.n))))
        nn = 0
        firstTerm = self.dlam_du(u)*(ivp(nn,self.l)-ivp(nn,self.l,2))*(self.zeta(nn)*self.Z(nn)+self.zeta(-nn)*self.Z(-nn))\
    +(iv(nn,self.l)-ivp(nn,self.l,1))*(self.dzetan_du(u,nn)*(self.Z(nn)+self.zeta(nn)*self.Zp(nn)) + self.dzetan_du(u,-nn)*(self.Z(-nn)+self.zeta(-nn)*self.Zp(-nn)))
        

        return self.pdByz_du(u)*self.Cyz() - self.Byz()*(sum + 0.5*firstTerm)




    def dKzz_du(self,u):
    

    
        sum= np.sum(self.dlam_du(u) * ivp(self.n,self.l)*( self.zeta(self.n)*self.Zp(self.n)+self.zeta(-self.n)*self.Zp(-self.n))\
                + iv(self.n,self.l)*( self.dzetan_du(u,self.n)*(self.Zp(self.n)+self.zeta(self.n)*self.Zpp(self.n)) + self.dzetan_du(u,-self.n)*( self.Zp(-self.n)+self.zeta(-self.n)*self.Zpp(-self.n) )  ))
    
        nn = 0
        firstTerm = self.dlam_du(u) * ivp(nn,self.l)*( self.zeta(nn)*self.Zp(nn)+self.zeta(-nn)*self.Zp(-nn))\
                + iv(nn,self.l)*( self.dzetan_du(u,nn)*(self.Zp(nn)+self.zeta(nn)*self.Zpp(nn)) + self.dzetan_du(u,-nn)*( self.Zp(-nn)+self.zeta(-nn)*self.Zpp(-nn) ) )
        
        
        return self.pdBzz_du(u)*(self.Kzz()-1) - self.Bzz()*(sum + 0.5*firstTerm)


    
    def pdBxx_du(self,u):
        return  1/self.X * H(u,"X") + 1/self.zeta(0)*self.dzetan_du(u,0) - self.dlam_du(u) - 1/self.l*self.dlam_du(u)

    def pdBxy_du(self,u):
        return  1/self.X * H(u,"X") + 1/self.zeta(0)*self.dzetan_du(u,0) - self.dlam_du(u) 

    def pdBxz_du(self,u):
        return  1/self.X * H(u,"X") - self.dlam_du(u) - 1/self.Yabs*H(u,"Y") - 1/self.Nz*H(u,"Nz") - 1/self.l*self.dlam_du(u)

    def pdByy_du(self,u):
        return  self.pdBxx_du(u)
    
    def pdByz_du(self,u):
        return  1/self.X * H(u,"X") - self.dlam_du(u) - 1/self.Yabs*H(u,"Y") - 1/self.Nz*H(u,"Nz") 

    def pdBzz_du(self,u):
        return  self.pdBxy_du(u)








