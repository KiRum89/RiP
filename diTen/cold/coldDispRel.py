import numpy as np
def sinTh(Nz,Nx):
    #
    return Nx/np.sqrt( Nx**2 + Nz**2 )

def cosTh(Nz,Nx):
    
    return Nz/np.sqrt( Nx**2 + Nz**2 )


def G(X,Yabs,Nz,Nx):
    return np.sqrt( Yabs**4*sinTh(Nz,Nx)**4 + 4*(1-X)**2*Yabs**2*cosTh(Nz,Nx)**2 )



def Q(X,Yabs,Nz,Nx,m):



    if m=="Om":
        return 2*(1-X)-Yabs**2*sinTh(Nz,Nx)**2+G(X,Yabs,Nz,Nx)
    if m=="Xm" :
        return 2*(1-X)-Yabs**2*sinTh(Nz,Nx)**2-G(X,Yabs,Nz,Nx)



def cold_disp_rel(X,Yabs,Nz,Nx,m):
    return 1-Nx**2-Nz**2 - 2*X*(1-X)/Q(X,Yabs,Nz,Nx,m)

def cold_disp_rel_Nx2(X,Yabs,Nz,Nx,m):
    return 1-Nz**2 - 2*X*(1-X)/Q(X,Yabs,Nz,Nx,m)


def coldDispRel_pol(X,Yabs,Nz,Nx):
    A = S(X,Yabs)
    B =  -( P(X,Yabs)*S(X,Yabs) + R(X,Yabs)*L(X,Yabs) ) + ( P(X,Yabs) + S(X,Yabs) )*Nz**2
    C = P(X,Yabs)*Nz**4 + P(X,Yabs)*R(X,Yabs)*L(X,Yabs) - ( 2*P(X,Yabs)*S(X,Yabs) )*Nz**2
    #print "A,B,C - cold",A,B,C
    return A*Nx**4 + B*Nx**2 + C

def coldDispRel_pol2(X,Yabs,Nz,Nx):
    cos2th= Nz**2/(Nx**2+Nz**2)
    sin2th= Nx**2/(Nx**2+Nz**2)
    A = S(X,Yabs)*sin2th + P(X,Yabs)*cos2th
    B = R(X,Yabs)*L(X,Yabs)*sin2th + P(X,Yabs)*S(X,Yabs)*(1+cos2th)
    C = P(X,Yabs)*R(X,Yabs)*L(X,Yabs)
    
    return A*( Nx**2+Nz**2 )**2 - B*( Nx**2+Nz**2 ) + C



def S(X,Yabs):
    return 1 - X/(1-Yabs**2)
def D(X,Yabs):
    return -Yabs*X/(1-Yabs**2)
def P(X,Yabs):
    return 1-X
def R(X,Yabs):
    return 1-X/(1-Yabs)
def L(X,Yabs):
    return 1-X/(1+Yabs)




def coldNx_Xm(X,Yabs,Nz):
    A = S(X,Yabs)
    B =  -( P(X,Yabs)*S(X,Yabs) + R(X,Yabs)*L(X,Yabs) ) + ( P(X,Yabs) + S(X,Yabs) )*Nz**2
    C = P(X,Yabs)*Nz**4 + P(X,Yabs)*R(X,Yabs)*L(X,Yabs) - ( 2*P(X,Yabs)*S(X,Yabs) )*Nz**2
    #print A,B,C
    return np.sqrt((-B-np.sqrt(B**2-4*A*C))/2/A)

#turning point when + and - modes concide
def turningPoint(X,Yabs,Nz):
    A = S(X,Yabs)
    B =  -( P(X,Yabs)*S(X,Yabs) + R(X,Yabs)*L(X,Yabs) ) + ( P(X,Yabs) + S(X,Yabs) )*Nz**2
    C = P(X,Yabs)*Nz**4 + P(X,Yabs)*R(X,Yabs)*L(X,Yabs) - ( 2*P(X,Yabs)*S(X,Yabs) )*Nz**2
    return B**2 - 4*A*C


def C(X,Yabs,Nz):
    C = P(X,Yabs)*Nz**4 + P(X,Yabs)*R(X,Yabs)*L(X,Yabs) - ( 2*P(X,Yabs)*S(X,Yabs) )*Nz**2
    return C

def coldNx_Om(X,Yabs,Nz):
    A = S(X,Yabs)
    B =  -( P(X,Yabs)*S(X,Yabs) + R(X,Yabs)*L(X,Yabs) ) + ( P(X,Yabs) + S(X,Yabs) )*Nz**2
    C = P(X,Yabs)*Nz**4 + P(X,Yabs)*R(X,Yabs)*L(X,Yabs) - ( 2*P(X,Yabs)*S(X,Yabs) )*Nz**2
    #print A,B,C
    return np.sqrt((-B+np.sqrt(B**2-4*A*C))/2/A)

