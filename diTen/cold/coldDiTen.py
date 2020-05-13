class ColdDiTen:

    def __init__(self,X,Yabs):
        self.X = X
        self.Yabs = Yabs

    def S(self):
        return 1 - self.X/(1-self.Yabs**2)

    def D(self):
        return -self.Yabs*self.X/(1-self.Yabs**2)

    def P(self):
        return 1-self.X

    def R(self):
        return 1-self.X/(1-self.Yabs)

    def L(self):
        return 1-self.X/(1+self.Yabs)
