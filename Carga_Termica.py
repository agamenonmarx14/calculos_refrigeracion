import math
from CoolProp.HumidAirProp import HAPropsSI
from Propiedades_Alimentos import PropAlimentos


class Cargas:
    def factor_usoDePuerta(P, theta_p, theta_o, theta_d):
        '''P= numero de transitos por puerta; theta_p= tiempo que tarda en abrir/cerrar puerta(s); theta_o= tiempo en que permanece abierta(min); theta_d= Intervalo de tiempo considerado en evaluación (h)'''
        Dt = (P*theta_p + 60*theta_o)/(3600*theta_d)
        return Dt

    def infiltracion(A, Ti, Tr, rhi, rhr, H, Dt, E=0.0):
        '''A= área de puerta; Ti= temperatura en espacio adyacente(°C); Tr= temperatura espacio refrigerado(°C); rhi= humedad relativa en espacio adyacente(0 a 1); rhr= humedad relativa en espacio refrigerado(0 a 1); H= altura de puerta(m); Dt= factor de tiempo de abertura de puerta(0 a 1); E= Efectividad de protección de puerta (0 a 1);\nRetorna carga térmica por infiltracion en Kw'''
        TD=Ti - Tr
        if TD <= 11:
            Df = 1.1
        else:
            Df = 0.8
        Ti = Ti + 273.15
        Tr = Tr + 273.15
        hi = HAPropsSI('H','T', Ti,'P',101325,'R', rhi)/1000
        hr = HAPropsSI('H','T', Tr,'P',101325,'R', rhr)/1000
        rho_i = 1/HAPropsSI('V','T', Ti,'P',101325,'R', rhi)
        rho_r = 1/HAPropsSI('V','T', Tr,'P',101325,'R', rhr)
        Fm = (2/(1+((rho_r/rho_i)**(1/3))))**1.5
        A = (1-(rho_i/rho_r))**0.5
        B = (9.81*H)**0.5
        q = 0.221*A*(hi - hr)*rho_r*A*B*Fm
        qt = q*Dt*Df*(1-E)
        return round(qt, 2)

    def carga_producto(producto, Tinicial, Tfinal, masa, t=1):
        '''producto: Producto a evaluar; Tinicial: Temperatura inicial(°C);\nTfinal: Temperatura final(°C); masa: masa de producto(Kg); t: Periodo evaluado(horas);\nRetorna carga térmica en Kw'''
        Ts = PropAlimentos.temperatura_inicioCongelacion(producto)
        c1 = PropAlimentos.Cp_SobreCong(producto, Tinicial)
        hs = PropAlimentos.Entalpia_Cong(producto, Ts)
        hf = PropAlimentos.Entalpia_BajoCong(producto, Tfinal)
        Q1 = m*c1*(Tinicial - Ts)
        Q2 = m*hs
        Q3 = m*(hs - hf)
        Q_producto = (Q1 + Q2 + Q3)/(3600*t)
        return round(Q_producto, 2)

    def carga_ocupacion(Tcamara, n):
        '''Tcamara: Temperatura de espacio refrigerado(°C); n: Número de personas.'''
        q = n*(272 - 6*Tcamara)
        return Q_ocupacion
