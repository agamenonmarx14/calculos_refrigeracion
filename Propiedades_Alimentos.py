import pandas as pd
import math

df = pd.read_csv('Propiedades.csv', index_col=0)

class PropAlimentos:
    def temperatura_inicioCongelacion(producto):
        t_congelacion = df['Ti'][producto]
        return t_congelacion

    def fraccion_agua(producto):
        return (df['xwo'][producto])/100

    def Tinicial_Cong(producto):
        return df['Ti'][producto]

    def agua_ligada(producto):
        xp = df['xp'][producto]
        xb = 0.4*xp
        return xb
    
    def fraccion_hielo_Tchigeov(producto, t):
        tf = df['Ti'][producto]
        xwo = df['xwo'][producto]/100
        if t <= tf:
            xice = 1.105*xwo/(1 + 0.7138/math.log(tf - t +1))
        else:
            xice = 0
        xu = xwo - xice
        return{'xice': xice, 'xu': xu}
    
    def fraccion_hielo_Miles(producto, t):
        tf = df['Ti'][producto]
        xwo = df['xwo'][producto]/100
        xb = 0.4*df['xp'][producto]/100
        if t <= tf:
            xice = (xwo-xb)*(1 - (tf/t))
        else:
            xice = 0
        xu = xwo - xice
        return{'xice': xice, 'xu': xu}
    
    def Cp_SobreCong(producto, t):
        if df['Ti'][producto] < t:
            xwo = df['xwo'][producto]/100
            xp = df['xp'][producto]/100
            xf = df['xf'][producto]/100
            xa = df['xa'][producto]/100
            xfb = df['xfb'][producto]/100
            xc = df['xc'][producto]/100

            cp = 2.0082 + 1.2089e-3 *t - 1.3129e-6* t**2
            cf = 1.9842 + 1.4733e-3 *t - 4.8008e-6 * t**2
            ca = 1.0926 + 1.8896e-3 *t - 3.6817e-6 *t**2
            cfb = 1.8459 + 1.8306e-3 * t - 4.6509e-6 *t**2
            cw = 4.1762 - 9.0864e-5 *t + 5.4731e-6 *t**2
            cc = 2.0141e-1 + 1.3874e-3*t - 4.3312e-6*t**2
            c = cw * xwo + cp * xp + cf*xf + ca*xa + cfb*xfb + cc*xc
        else:
            c = 0
        return c
    
    def Cp_BajoConge(producto, t):
        tf = df['Ti'][producto]
        xp = df['xp'][producto]/100
        xwo = df['xwo'][producto]/100
        L0 = df['Latent'][producto]
        ca = 1.55 + 1.26*(1-xwo) - (xwo - 0.4*xp)*L0*tf/t**2
        return ca
    
    def Entalpia_SobreCong(producto, t):
        tf = df['Ti'][producto]
        tr = -40
        xwo = df['xwo'][producto]/100
        xs = 1 - xwo
        xb = 0.4*df['xp'][producto]/100
        L0 = df['Latent'][producto]
        Hf = (tf - tr)*(1.55 + 1.26*xs - ((xwo - xb)*L0*tf/(tr*tf)))
        H = Hf + (t - tf)*(4.19 - 2.3*xs - 0.628*xs**3)
        return H
    
    def Entalpia_Cong(producto):
        tf = df['Ti'][producto]
        tr = -40
        xwo = df['xwo'][producto]/100
        xs = 1 - xwo
        xb = 0.4*df['xp'][producto]/100
        L0 = df['Latent'][producto]
        H = (tf - tr)* (1.55 + 1.26*xs - ( ((xwo - xb)*L0*tf)/(tf*tr) ) )
        return H
    
    def Entalpia_BajoCong(producto, t):
        tr = -40
        tf = df['Ti'][producto]
        xwo = df['xwo'][producto]/100
        xs = 1 - xwo
        xb = 0.4*df['xp'][producto]/100
        L0 = df['Latent'][producto]
        H = (t - tr)*(1.55 + 1.26*xs - ((xwo - xb)*L0*tf/(tr*t)))
        return H
    
    def Conductividad_Paralela(producto, t):
        xp=df['xp'][producto]/100
        xf=df['xf'][producto]/100
        xfb=df['xfb'][producto]/100
        xa=df['xa'][producto]/100
        xc=df['xc'][producto]/100
        xice= PropAlimentos.fraccion_hielo_Miles(producto, t)['xice']
        xwo=PropAlimentos.fraccion_hielo_Miles(producto, t)['xu']
        kp=1.7881e-1 + (1.1958e-3) * t - (2.7178e-6) * t**2
        kf=1.8071e-1 - (2.7604e-4) * t - (1.7749e-7) * t**2
        kc=2.0141e-1 + (1.3874e-3) * t - (4.3312e-6) * t**2
        kfb=1.8331e-1 + (1.2497e-3) * t - (3.1683e-6) * t**2
        ka=3.2962e-1 + (1.4011e-3) * t - (2.9069e-6) * t**2
        kice = 2.2196 - (6.2489e-3) * t + (1.0154e-4) * t**2
        kwo=5.7109e-1 + (1.7625e-3) * t - (6.7036e-6) * t**2
        Dp = 1.3299e3 - (5.1840e-1)*t
        Df = 9.2559e2 - (4.1757e-1)*t
        Da = 2.4238e3 - (2.8063e-1)*t
        Dc = 1.5991e3 - (3.1046e-1)*t
        Dfb = 1.3115e3 - (3.6589e-1)*t
        Dice = 9.1689e2 - (1.3071e-1)*t
        Dw = 9.9718e2 + (3.1439e-3)*t - (3.7574e-3)*t**2
        x_D = ((xp/Dp) + (xf/Df) + (xfb/Dfb) + (xa/Da) + (xc/Dc) + 
               (xice/Dice) + (xwo/Dw))
        k = ((kp*(xp/Dp)/x_D) + (kf*(xf/Df)/x_D) + (kfb*(xfb/Dfb)/x_D) +
             (kc*(xc/Dc)/x_D) + (ka*(xa/Da)/x_D) + (kwo*(xwo/Dw)/x_D)+
             (kice*(xice/Dice)/x_D))
        return k

    def Conductividad_Perpendicular(producto, t):
        xp=df['xp'][producto]/100
        xf=df['xf'][producto]/100
        xfb=df['xfb'][producto]/100
        xa=df['xa'][producto]/100
        xc=df['xc'][producto]/100
        xice= PropAlimentos.fraccion_hielo_Tchigeov(producto, t)['xice']
        xwo=PropAlimentos.fraccion_hielo_Tchigeov(producto, t)['xu']
        kp=1.7881e-1 + (1.1958e-3) * t - (2.7178e-6) * t**2
        kf=1.8071e-1 - (2.7604e-3) * t - (1.7749e-7) * t**2
        kc=2.0141e-1 + (1.3874e-3) * t - (4.3312e-6) * t**2
        kfb=1.8331e-1 + (1.2497e-3) * t - (3.1683e-6) * t**2
        ka=3.2962e-1 + (1.4011e-3) * t - (2.9069e-6) * t**2
        kice = 2.2196 - (6.2489e-3) * t + (1.0154e-4) * t**2
        kwo=5.7109e-1 + (1.7625e-3) * t - (6.7036e-6) * t**2
        Dp = 1.3299e3 - (5.1840e-1)*t
        Df = 9.2559e2 - (4.1757e-1)*t
        Da = 2.4238e3 - (2.8063e-1)*t
        Dc = 1.5991e3 - (3.1046e-1)*t
        Dfb = 1.3115e3 - (3.6589e-1)*t
        Dice = 9.1689e2 - (1.3071e-1)*t
        Dw = 9.9718e2 + (3.1439e-3)*t - (3.7574e-3)*t**2
        x_D = ((xp/Dp) + (xf/Df) + (xfb/Dfb) + (xa/Da) + (xc/Dc) + 
               (xice/Dice) + (xwo/Dw))
        k = 1 / (((xp/Dp)/x_D)/kp + ((xf/Df)/x_D)/kf + ((xfb/Dfb)/x_D)/kfb + ((xa/Da)/x_D)/ka + ((xc/Dc)/x_D)/kc + 
               ((xice/Dice)/x_D)/kice + ((xwo/Dw)/x_D)/kwo)
        return k
    
    def Densidad(producto, t):
        tf=df['Ti'][producto]
        xp=df['xp'][producto]/100
        xf=df['xf'][producto]/100
        xfb=df['xfb'][producto]/100
        xa=df['xa'][producto]/100
        xc=df['xc'][producto]/100
        if t==tf:
            xice= PropAlimentos.fraccion_hielo_Miles(producto, t)['xice']
            xwo=PropAlimentos.fraccion_hielo_Miles(producto, t)['xu']
        else:
            xice= PropAlimentos.fraccion_hielo_Tchigeov(producto, t)['xice']
            xwo=PropAlimentos.fraccion_hielo_Tchigeov(producto, t)['xu']
        #elif t<tf:
        Dp = 1.3299e3 - (5.1840e-1)*t
        Df = 9.2559e2 - (4.1757e-1)*t
        Da = 2.4238e3 - (2.8063e-1)*t
        Dc = 1.5991e3 - (3.1046e-1)*t
        Dfb = 1.3115e3 - (3.6589e-1)*t
        Dice = 9.1689e2 - (1.3071e-1)*t
        Dw = 9.9718e2 + (3.1439e-3)*t - (3.7574e-3)*t**2
        D = 1/(xp/Dp + xf/Df + xfb/Dfb + xa/Da + xc/Dc + xice/Dice + xwo/Dw)
        return D

    def Difusividad(producto, T):
        Tf=df['Ti'][producto]
        k = PropAlimentos.PropAlimentos.Conductividad_Paralela(producto, T)
        rho = PropAlimentos.Densidad(producto, T)
        if T > Tf:
            c=PropAlimentos.Cp_SobreCong(producto, T)
        else:
            c=PropAlimentos.Cp_BajoConge(producto, Ti)
        alpha = k/(rho*c)
        return round(2, alpha)
    
    def Omega(f, a, Bi):
        tol=0.001
        while True:
            u = abs(f(a) - 1)
            v = abs(Bi)
            if abs(u-v) <= tol:
                break
            else:
                a = a + 0.00001
        return(round(a, 3)) 

    def Tiempo_Enfriamiento(producto, Tm, Ti, T, dimensiones, forma, h):
        t = (Ti - T)/2
        rho = PropAlimentos.Densidad(producto, t)
        k = PropAlimentos.Conductividad_Perpendicular(producto, t)
        c = PropAlimentos.Cp_SobreCong(producto, T)*1000
        L =min(dimensiones)
        dimensiones.remove(L)
        Lbiot = L/2
        beta1 = min(dimensiones)/L
        beta2 = max(dimensiones)/L
        GParams = {
            'Losa Infinita':{'N': 1, 'p1':0, 'p2':0, 'p3':0, 'y1':1e20, 'y2':1e20, 'l':1},
            'Varilla Rectangular Infinita': {'N': 2, 'p1':0.75, 'p2':0, 'p3':-1, 'y1':(4*beta1)/math.pi, 'y2':1e20, 'l':(4*beta1)/math.pi},
            'Bloque': {'N': 3, 'p1':0.75, 'p2':0.75, 'p3':-1, 'y1':(4*beta1)/math.pi, 'y2':1.5*beta2, 'l':(4*beta1)/math.pi},
            'Cilindro Infinito': {'N': 2, 'p1':1.01, 'p2':0, 'p3':0, 'y1':1, 'y2':1e20, 'l':1},
            'Elipse Infinito': {'N': 2, 'p1':1.01, 'p2':0, 'p3':1, 'y1':beta1, 'y2':1e20, 'l':beta1},
            'Cilindro Achatado': {'N': 3, 'p1':1.01, 'p2':0.75, 'p3':-1, 'y1':1.225*beta1, 'y2':1.225*beta2, 'l':1.225*beta1},
            'Cilindro Corto': {'N': 3, 'p1':1.01, 'p2':0.75, 'p3':-1, 'y1':beta1, 'y2':1.5*beta2, 'l':beta1},
            'Esfera' : {'N': 3, 'p1':1.01, 'p2':1.24, 'p3':0, 'y1': 1, 'y2':1, 'l':1},
            'Elipsoide': {'N': 3, 'p1':1.01, 'p2':1.24, 'p3':1, 'y1':beta1, 'y2':beta2, 'l':beta1}
        }
        Bi = (h*Lbiot)/k
        E0 = 1.5 * ((beta1 + beta2 + (1+beta2)*beta1**2 + (1+beta1)*beta2**2)/(beta1*beta2*(beta1 + beta2 +1))) - (((beta1 - beta2)**2)**0.4)/15
        p1 = GParams[forma]['p1']
        p2 = GParams[forma]['p2']
        p3 = GParams[forma]['p3']
        Fbeta1 = (1/beta1**2) + 0.01*p3*math.exp(beta1 -(beta1**2)/6)
        Fbeta2 = (1/beta2**2) + 0.01*p3*math.exp(beta2 -(beta2**2)/6)
        Ei = 0.75 + p1*Fbeta1 + p2*Fbeta2
        E = (1.85 + (Bi**(4/3))) / ((Bi**(4/3)/Ei)+(1.85/E0))
        y1 = GParams[forma]['y1']
        y2 = GParams[forma]['y2']
        lamb = GParams[forma]['l']
        N = GParams[forma]['N']
        ji = 1.271 + 0.305*math.exp(0.172*y1 - 0.115*y1**2) + 0.425*math.exp(0.09*y2 - 0.128*y2**2)
        jc = ((Bi**1.35) + (1/lamb)) / (((Bi**1.35)/ji) + (1/lamb))
        jm = jc * (((1.5 + 0.69*Bi)/(1.5 + Bi))**N)
        Y = (Tm - T)/(Tm -Ti)
        f = lambda x: x/math.tan(x)
        omega=PropAlimentos.Omega(f, 0.001, Bi)
        theta = ((3*rho*c*(Lbiot**2))/(k*E*omega**2)*(math.log(jm/Y)))
        return (round(theta, 2))
    
    def Tiempo_Congelacion_EHTD(producto, Ti, Te, dimensiones, Tm, h, forma):

        '''producto: Producto a evaluar; Ti: Temperatura inicial del centro del producto(°C); \ndimensiones: Lista con las dimensiones del producto en metros;\nTm: Temperatura del espacio de Tiempo_Enfriamiento(°C); h: Coeficiente de transferencia por convección (W/(m²K);\nforma: "Esfera", "Placa Infinita", "Cilindro Infinito", "Bloque")'''

        Tf=df['Ti'][producto]
        Di=PropAlimentos.Densidad(producto=producto, t=Ti)
        De=PropAlimentos.Densidad(producto=producto, t=Te)
        Df=PropAlimentos.Densidad(producto=producto, t=Tf)
        Dff=PropAlimentos.Densidad(producto=producto, t=-40)
        kff=PropAlimentos.Conductividad_Paralela(producto, -40)
        Cff = PropAlimentos.Cp_BajoConge(producto, -40)
        if Ti < Tf:
            Ci=PropAlimentos.Cp_BajoConge(producto, Ti)
        else:
            Ci=PropAlimentos.Cp_SobreCong(producto, Ti)
        He = PropAlimentos.Entalpia_BajoCong(producto, Te)
        Hf = PropAlimentos.Entalpia_Cong(producto)

        Hvf = Df*Hf - De*He
        DCi = Ci*Di
        DCff = Dff*Cff
        dmenor = min(dimensiones)
        D=dmenor
        dimensiones.remove(dmenor)
        beta1=min(dimensiones)/D
        beta2=max(dimensiones)/D
        Bi = (h*D)/kff
        Pk = (DCi*(Ti - Tf))/Hvf
        Ste=(DCff*(Tf - Tm)) / Hvf

        def PPlaca_Infinita(Pk, Ste, Bi, D):
            if 0 <= D and D <= 0.12:
                if 10 <= h and h <=500:
                    if Ti <= 40:
                        if -45 <= Tm and Tm <= 15:
                            P= 0.5072 + 0.2018*Pk + Ste*(0.3224*Pk + (0.0105/Bi) + 0.0681)
                        else:
                            raise ValueError('Tm debe cumplir -45<=Tm<=-15')
                    else:
                        raise ValueError('Ti debe ser menor a 40 °C')
                else:
                    raise ValueError('h Debe cumplir 10<=h<=500')
            else:
                raise ValueError('D dede cumplir 0<=D<=0.12')
            return P
        
        def RPlaca_Infinita(Pk, Ste, Bi, D):
            if 0 <= D and D <= 0.12:
                R = 0.1684 + Ste*(0.274*Pk - 0.0135)
            else:
                raise ValueError('D dede cumplir 0<=D<=0.12')
            return R
        
        def PCilindro_Infinito(Pk, Ste, Bi):
            if 0.155 <= Ste and Ste <= 0.345:
                if 0.5 <= Bi and Bi <= 4.5:
                    if 0 <= Pk and Pk <= 0.55:
                        P = 0.3751+0.09999*Pk +Ste*(0.4008*Pk + (0.071/Bi)-0.5865)
                    else:
                        raise ValueError('Valor de Pk debe cumplir 0<=Pk<=0.55')
                else:
                    raise ValueError('Valor de Bi debe cumplir 0.5<=Bi<=4.5')
            else:
                raise ValueError('Valor Ste debe cumplir 0.155<=Ste<= 0.345')
            return P

        def RCilindro_Infinito(Pk, Ste, Bi):
            if 0.155 <= Ste and Ste <= 0.345:
                if 0.5 <= Bi and Bi <= 4.5:
                    if 0 <= Pk and Pk <= 0.55:
                        R = 0.0133 + Ste*(0.0415*Pk + 0.3957)
                    else:
                        raise ValueError('Valor de Pk debe cumplir 0<=Pk<=0.55')
                else:
                    raise ValueError('Valor de Bi debe cumplir 0.5<=Bi<=4.5')
            else:
                raise ValueError('Valor Ste debe cumplir 0.155<=Ste<= 0.345')  
            return R

        def PEsfera(Pk, Ste, Bi):
            if 0.155 <= Ste and Ste <= 0.345:
                if 0.5 <= Bi and Bi <= 4.5:
                    if 0 <= Pk and Pk <= 0.55:
                        P = 0.1084 + 0.0924*Pk + Ste*(0.231 - (0.3114/Bi) + 0.6739)
                    else:
                        raise ValueError('Valor de Pk debe cumplir 0<=Pk<=0.55')
                else:
                    raise ValueError('Valor de Bi debe cumplir 0.5<=Bi<=4.5')
            else:
                raise ValueError('Valor Ste debe cumplir 0.155<=Ste<= 0.345')
            return P

        def REsfera(Pk, Ste, Bi):
            if 0.155 <= Ste and Ste <= 0.345:
                if 0.5 <= Bi and Bi <= 4.5:
                    if 0 <= Pk and Pk <= 0.55:
                        R = 0.0784 + Ste*(0.0386*Pk - 0.1694)
                    else:
                        raise ValueError('Valor de Pk debe cumplir 0<=Pk<=0.55')
                else:
                    raise ValueError('Valor de Bi debe cumplir 0.5<=Bi<=4.5')
            else:
                raise ValueError('Valor Ste debe cumplir 0.155<=Ste<= 0.345')
            return R                

        def PBloque(Pk, Ste, Bi, beta1, beta2):
            if 0.155 <= Ste and Ste <= 0.345:
                if 0.5 <= Bi and Bi <= 4.5:
                    if 0 <= Pk and Pk <= 0.55:
                        if 1 <= beta1 and beta1 <= 4:
                            if 1 <= beta1 and beta1 <= 4:
                                P1 = (beta1*beta2)/(2*(beta1*beta2 + beta1 + beta2))
                                P2 = P1*(1.026 + 0.5808*Pk + Ste*(0.2296*Pk + (0.0182/Bi) + 0.105))
                                P = P2 + P1*(0.1136 + Ste*(5.766*P1 - 1.242))
                            else:
                                raise ValueError('Valor de beta2 debe cumplir 1<=beta2<=4')
                        else:
                            raise ValueError('Valor de beta1 debe cumplir 1<=beta1<=4')
                    else:
                        raise ValueError('Valor de Pk debe cumplir 0<=Pk<=0.55')
                else:
                    raise ValueError('Valor de Bi debe cumplir 0.5<=Bi<=4.5')
            else:
                raise ValueError('Valor Ste debe cumplir 0.155<=Ste<= 0.345')
            return P

        def RBloque(Pk, Ste, Bi, beta1, beta2):
            if 0.155 <= Ste and Ste <= 0.345:
                if 0.5 <= Bi and Bi <= 4.5:
                    if 0 <= Pk and Pk <= 0.55:
                        if 1 <= beta1 and beta1 <= 4:
                            if 1 <= beta1 and beta1 <= 4:
                                Q=4*((beta1-beta2)*(beta1-1) +((beta2-1)**2))**0.5
                                r = (1/3)*(beta1 +  beta2 + 1 + ((beta1-beta2)*(beta1-1)+(beta2-1)**2)**0.5)
                                s = (1/3)*(beta1 + beta2 + 1 - ((beta1-beta2)*(beta1-1)+(beta2-1)**2)**0.5)
                                e1 = (r - 1)*(beta1 - r)*(beta2 - r)*math.log((r/(r - 1)))
                                e2 = (s - 1)*(beta1 - s)*(beta2 - s)*math.log((s/(s-1)))
                                R1 = ((1/(2*Q))*(e1 - e2)) + (1/72)*(2*beta1 + 2*beta2 -1)
                                R2 = R1*(1.202 + Ste*(3.41*Pk + 0.7336))
                                R = R2 + R1*(0.7334 + Ste*(49.86*R1 - 2.9))
                            else:
                                raise ValueError('Valor de beta2 debe cumplir 1<=beta2<=4')
                        else:
                            raise ValueError('Valor de beta1 debe cumplir 1<=beta1<=4')
                    else:
                        raise ValueError('Valor de Pk debe cumplir 0<=Pk<=0.55')
                else:
                    raise ValueError('Valor de Bi debe cumplir 0.5<=Bi<=4.5')
            else:
                raise ValueError('Valor Ste debe cumplir 0.155<=Ste<= 0.345')
            return R

        RDict = {'Esfera': REsfera, 'Placa Infinita': RPlaca_Infinita, 'Cilindro Infinito': RCilindro_Infinito, 'Bloque': RBloque}    
        PDict = {'Esfera': PEsfera, 'Placa Infinita': PPlaca_Infinita, 'Cilindro Infinito': PCilindro_Infinito, 'Bloque': PBloque}

        if forma == 1 or forma == 'Bloque':
            R = RDict['Bloque'](Pk, Ste, Bi, beta1, beta2)
            P = PDict['Bloque'](Pk, Ste, Bi, beta1, beta2)
        elif forma == 2 or forma == 'Esfera':
            R = RDict['Esfera'](Pk, Ste, Bi)
            P = PDict['Esfera'](Pk, Ste, Bi)
        elif forma == 3 or forma == 'Cilindro Infinito':
            R = RDict['Cilindro Infinito'](Pk, Ste, Bi)
            P = PDict['Cilindro Infinito'](Pk, Ste, Bi)
        elif forma == 4 or forma == 'Placa Infinita':
            R = RDict['Cilindro Infinito'](Pk, Ste, Bi, h, D)
            P = PDict['Cilindro Infinito'](Pk, Ste, Bi, h, D)
        
        theta = (((1000*Hvf/(Tf - Tm)) * ( (P*D/h) + (R*D**2)/kff )) * (1 - (1.65*Ste/kff)*math.log((Te - Tm)/(-10 - Tm))))
        return round(theta, 2)
