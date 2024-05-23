import CoolProp.CoolProp as CP
from CoolProp.HumidAirProp import HAPropsSI
import math

class ProcesosAire:
    def __init__(self, Tdb, HR=100, P=101.325):
        self.P = P*1000
        self.Tdb = Tdb + 273.15
        self.HR = HR/100
        self.densidad = CP.PropsSI('D', 'T', self.Tdb, 'P', self.P, 'Air')

    def enfriamiento_mixto(self, Tdbi, HRi, Tdbf, HRf):
        Tdbi = Tdbi + 273.15
        Tdbf = Tdbf + 273.15
        HRi = HRi/100
        HRf = HRf/100
        Hi = HAPropsSI('H','T', Tdbi,'P',101325,'R', HRi)
        Hf = HAPropsSI('H','T', Tdbf,'P',101325,'R', HRf)
        deltaH = (Hf - Hi)/1000
        return deltaH

class CompresorReciprocante:
    def __init__(self, ref, vol=1.0, ni=1.0, nv=1.0):
        '''
        ref = refrigerante. Selecciona el refrigerante o gas a analizar. Opciones:/n
        R744, R717, Aire, Helio, Nitrógeno, R12, R1234yf, R134a, R22, R32, R404A, R410A./n
        vol = Desplazamiento volumétrico en m³/h. Default: 1.0 m³/h./n
        '''
        self.ref = ref
        self.vol = vol
        self.refrigerantes = {
            'R744':'CarbonDioxide', 'R717':'Ammonia', 'Aire':'Air', 'Helio':'Helium',
            'Nitrógeno':'Nitrogen', 'Nitrogeno':'Nitrogen', 'R12':'R12', 'R1234yf':'R1234yf',
            'R134a':'R134a', 'R22':'R22', 'R32':'R32', 'R404A':'R404A',
            'R410A':'R410A'
            }

    def eficiencia_volumetrica(self, claro, recalentamiento_aspiracion, ps, pd, pld=0.01, pls=0.01, recalentamiento_camara=0.0):
        '''
        ps = Presión de succión (Bar)./n
        pd = Presión de descarga (Bar)/N
        claro = Porcentaje de claro de válvulas dividido entre cien./n
        pls = Caída de presión en válvulas de descarga(Bar)./n
        pld = Caída de presión en válvulas de aspiración(Bar)./n
        '''
        refrigerante = self.refrigerantes[self.ref]
        Pb = (ps+1.01325)*100000
        Pc = (pd+1.01325)*100000
        ts = CP.PropsSI('T', 'P', Pb, 'Q', 0, ) + recalentamiento_aspiracion
        d1 = CP.PropsSI('D', 'T|gas', Pb, 'P', ps, refrigerante)
        cp = CP.PropsSI('Cpmass', 'T|gas', ts, 'P', Pb, refrigerante)
        cv = CP.PropsSI('Cvmass', 'T|gas', ts, 'P', Pb, refrigerante)
        Pb = (ps - pls +1.01325)*100000
        Pc = (pd + pld +1.01325)*100000
        d2 = CP.PropsSI('D', 'T|gas', (ts+recalentamiento_camara), 'P', Pb, refrigerante)
        n = cp/cv
        nv = (1 + claro - claro*((Pc/Pb)**(1/n)))*(d2/d1)
        return nv

    def flujo_masico(self, ps, claro, ssh, nv=1.0, i_sh=10):
        try:
            ps = (ps+1.01325)*100000
            ts = CP.PropsSI('T', 'P', ps, 'Q', 0, self.refrigerantes[self.ref]) + ssh
            densidad_vapor = CP.PropsSI('D', 'T|gas', ts, 'P', ps, self.refrigerantes[self.ref])
            flujo = (self.vol*nv*densidad_vapor)/3600
        except:
            print('Datos proporcionados son erróneos')
        return flujo

    def trabajo_compReciprocante(self, ps, pd, recalentamiento_aspiracion, flujo_refrigerante, recalentamiento_camara=0.0, ni=1.0):
        refrigerante = self.refrigerantes[self.ref]
        Pb = (ps+1.01325)*100000
        Pc = (pd+1.01325)*100000
        ts = CP.PropsSI('T', 'P', Pb, 'Q', 1, refrigerante) + recalentamiento_aspiracion;
        t_camara = ts + recalentamiento_camara
        cp = CP.PropsSI('Cpmass', 'T|gas', ts, 'P', Pb, refrigerante)
        cv = CP.PropsSI('Cvmass', 'T|gas', ts, 'P', Pb, refrigerante)
        Pb = Pb - pls*100000
        Pb = Pc + pld*100000
        vb = 1/CP.PropsSI('D', 'T|gas', t_camara, 'P', Pb, refrigerante)
        n = cp/cv
        w = (n/(n-1)) * (Pb*vb) * (((Pc/Pb)**((n-1)/n)) -1)
        w = w*flujo_refrigerante/(ni*1000)
        return w

    def eficiencia_isentropica(self, ps, pd, td, desplazamiento):
        ps = (ps + 1.01325)*100000
        pd = (pd + 1.01325)*100000
        td = td + 273.15
        ts = CP.PropsSI('T', 'P', ps, 'Q', 0, self.refrigerantes[self.ref]) + ssh
        hi = CP.PropsSI('H','T', ts, 'P', ps, self.refrigerantes[self.ref])
        s = CP.PropsSI('S','T', ts,'P', ps, self.refrigerantes[self.ref])
        hs = CP.PropsSI('H','S', s,'P', pd, self.refrigerantes[self.ref])
        hr = CP.PropsSI('H','T', td,'P', pd, self.refrigerantes[self.ref])
        self.ni = (hs - hi)/(hr - hi)
        return self.ni

class CompresorScroll:
    def __init__(self, ref, vol=1.0, ni=1.0, nv=1.0):
        '''
        ref = refrigerante. Selecciona el refrigerante o gas a analizar. Opciones:/n
        R744, R717, Aire, Helio, Nitrógeno, R12, R1234yf, R134a, R22, R32, R404A, R410A./n
        desplazamiento = Desplazamiento volumétrico en m³/h. Default: 1.0 m³/h./n
        '''
        self.ref = ref
        self.desplazamiento = desplazamiento
        self.refrigerantes = {
            'R744':'CarbonDioxide', 'R717':'Ammonia', 'Aire':'Air', 'Helio':'Helium',
            'Nitrógeno':'Nitrogen', 'Nitrogeno':'Nitrogen', 'R12':'R12', 'R1234yf':'R1234yf',
            'R134a':'R134a', 'R22':'R22', 'R32':'R32', 'R404A':'R404A',
            'R410A':'R410A'
            }

    def trabajo_compresion(self, ps, pd, recalentamiento_aspiracion, desplazamiento, Vr, nv=1):
        '''
        ps: Presión de aspiración(bar); pd: Presión de condensación(bar); desplazamiento: Desplazamiento volumétrico(m³/h); Vr: Relación de volumen; nv: Eficiencia volumétrica(Default=1.0).\nResultado: Potencia consumida en Kw
        '''
        Pb = (ps+1.01325)*100000;#Presión de aspiración
        Pc = (pd+1.01325)*100000;#Presión de descarga
        ts = CP.PropsSI('T', 'P', Pb, 'Q', 1, self.refrigerantes[self.ref]) + recalentamiento_aspiracion;#Temperatura de succión.
        cp = CP.PropsSI('Cpmass', 'T|gas', ts, 'P', Pb, self.refrigerantes[self.ref])
        cv = CP.PropsSI('Cvmass', 'T|gas', ts, 'P', Pb, self.refrigerantes[self.ref])
        D = CP.PropsSI('D', 'T|gas', ts, 'P', Pb, self.refrigerantes[self.ref]);#Volumen específico succión
        vb = 1/D
        n = cp/cv
        w = (nv*(desplazamiento/3600)*D)*((n/(n-1)) * (Pb*vb)*((Vr**(n-1)) -1) - (vb/Vr)*(Pb*(Vr**n) - Pc))
        return (w/1000)

class EvaporadorDX:
        def __init__(self, ref, pc, pe, sh, sc):
            self.ref = ref
            self.pc = pc
            self.pe = pe
            self.sh = sh
            self.sc = sc
            self.refrigerantes = {
                'R744':'CarbonDioxide', 'R717':'Ammonia', 'Aire':'Air', 'Helio':'Helium',
                'Nitrógeno':'Nitrogen', 'Nitrogeno':'Nitrogen', 'R12':'R12', 'R1234yf':'R1234yf',
                'R134a':'R134a', 'R22':'R22', 'R32':'R32', 'R404A':'R404A',
                'R410A':'R410A'
                }

        def capacidad_frigorifica(self, flujo_refrigerante):
            pc = (self.pc+1.01325)*100000
            pe = (self.pe+1.01325)*100000
            tc = CP.PropsSI('T', 'P', pc, 'Q', 0, self.refrigerantes[self.ref])
            te = CP.PropsSI('T', 'P', pe, 'Q', 0, self.refrigerantes[self.ref])
            hi = CP.PropsSI('H','T|liquid', (tc-self.sc),'P', pc , self.refrigerantes[self.ref])
            hf = CP.PropsSI('H','T|gas', (te+self.sh), 'P', pe, self.refrigerantes[self.ref])
            deltah = (hf - hi)/1000
            capacidad = flujo_refrigerante*deltah
            return capacidad
