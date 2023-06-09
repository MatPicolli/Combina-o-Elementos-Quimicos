from os import system

class Elemento:
    def __init__(self, simbolo, nome, valencia):
        self.simbolo = simbolo
        self.nome = nome
        self.valencia = valencia

class Combinacao:
    def __init__(self, *elementos):
        # Ordena os elementos por valência e, em seguida, por ordem alfabética
        self.elementos = sorted(elementos, key=lambda e: e.valencia, reverse=True)

    def _calcula_prefixos(self, valencia1, valencia2):
        prefixos = {
            1: 'mono', 2: 'di', 3: 'tri', 4: 'tetra', 5: 'penta', 6: 'hexa', 7: 'hepta', 8: 'octa', 9: 'nona', 10: 'deca'
        }
        return (prefixos[abs(valencia1)], prefixos[abs(valencia2)])

    def nome(self):
        # Verifica se algum elemento tem valência 0 e se possui valências específicas
        valencia_zero = False
        for elem in self.elementos:
            if elem.valencia == 0:
                if not elem.valencias_especificas:
                    valencia_zero = True
                else:
                    for val, simb in elem.valencias_especificas:
                        if any(e.simbolo == simb for e in self.elementos):
                            elem.valencia = val
                            break

        if valencia_zero:
            return "Um ou mais elementos possuem valência 0 e não formam compostos nessa combinação."

        # Identifica o elemento com menor valência e o com maior valência
        elem_menor_valencia = min(self.elementos, key=lambda e: abs(e.valencia))
        elem_maior_valencia = max(self.elementos, key=lambda e: abs(e.valencia))

        # Troca a ordem dos elementos, se necessário
        if elem_menor_valencia != self.elementos[0]:
            self.elementos.reverse()

        # Calcula os prefixos dos nomes dos elementos
        prefixos = self._calcula_prefixos(self.elementos[0].valencia, self.elementos[1].valencia)

        # Monta o nome do composto
        nome_elem1 = self.elementos[0].nome
        nome_elem2 = self.elementos[1].nome

        if self.elementos[1].valencia < 0:
            nome_elem2 = nome_elem2[:-1] + 'eto'
        else:
            nome_elem2 = prefixos[1] + nome_elem2

        nome_composto = f"{nome_elem1}eto de {nome_elem2}"
        return nome_composto

    def formula(self):
        # Cria uma lista com os símbolos dos elementos da combinação
        simbolos = [elem.simbolo for elem in self.elementos]

        # Cria uma lista com as valências dos elementos da combinação
        valencias = [abs(elem.valencia) for elem in self.elementos]

         # Cria a fórmula química combinando os símbolos e valências dos elementos
        formula_quimica = ''.join([s + (str(v) if v > 1 else '') for s, v in zip(simbolos, valencias)])

        return formula_quimica


# cria os elementos químicos
H = Elemento('H', 'hidrogênio', 1)
He = Elemento('He', 'hélio', 0)
Li = Elemento('Li', 'lítio', 1)
Be = Elemento('Be', 'berílio', 2)
B = Elemento('B', 'boro', 3)
C = Elemento('C', 'carbono', 4)
N = Elemento('N', 'nitrogênio', -3)
O = Elemento('O', 'oxigênio', -2)
F = Elemento('F', 'flúor', -1)
Ne = Elemento('Ne', 'neônio', 0)
Na = Elemento('Na', 'sódio', 1)
Mg = Elemento('Mg', 'magnésio', 2)
Al = Elemento('Al', 'alumínio', 3)
Si = Elemento('Si', 'silício', 4)
P = Elemento('P', 'fósforo', -3)
S = Elemento('S', 'enxofre', -2)
Cl = Elemento('Cl', 'cloro', -1)
Ar = Elemento('Ar', 'argônio', 0)
K = Elemento('K', 'potássio', 1)
Ca = Elemento('Ca', 'cálcio', 2)
Sc = Elemento('Sc', 'escândio', 3)
Ti = Elemento('Ti', 'titânio', 4)
V = Elemento('V', 'vanádio', 5)
Cr = Elemento('Cr', 'cromo', 6)
Mn = Elemento('Mn', 'manganês', 7)
Fe = Elemento('Fe', 'ferro', 2)
Co = Elemento('Co', 'cobalto', 3)
Ni = Elemento('Ni', 'níquel', 2)
Cu = Elemento('Cu', 'cobre', 2)
Zn = Elemento('Zn', 'zinco', 2)
Ga = Elemento('Ga', 'gálio', 3)
Ge = Elemento('Ge', 'germânio', 4)
As = Elemento('As', 'arsênio', -3)
Se = Elemento('Se', 'selênio', -2)
Br = Elemento('Br', 'bromo', -1)
Kr = Elemento('Kr', 'criptônio', 0)
Rb = Elemento('Rb', 'rubídio', 1)
Sr = Elemento('Sr', 'estrôncio', 2)
Y = Elemento('Y', 'ítrio', 3)
Zr = Elemento('Zr', 'zircônio', 4)
Nb = Elemento('Nb', 'nióbio', 5)
Mo = Elemento('Mo', 'molibdênio', 6)
Tc = Elemento('Tc', 'tecnécio', 7)
Ru = Elemento('Ru', 'rutênio', 4)
Rh = Elemento('Rh', 'ródio', 3)
Pd = Elemento('Pd', 'paládio', 2)
Ag = Elemento('Ag', 'prata', 1)
Cd = Elemento('Cd', 'cádmio', 2) 
In = Elemento('In', 'índio', 3)
Sn = Elemento('Sn', 'estanho', 4)
Sb = Elemento('Sb', 'antimônio', -3)
Te = Elemento('Te', 'telúrio', -2)
I = Elemento('I', 'iodo', -1)
Xe = Elemento('Xe', 'xenônio', 0)
Cs = Elemento('Cs', 'césio', 1)
Ba = Elemento('Ba', 'bário', 2)
La = Elemento('La', 'lantânio', 3)
Ce = Elemento('Ce', 'cério', 4)
Pr = Elemento('Pr', 'praseodímio', 3)
Nd = Elemento('Nd', 'neodímio', 3)
Pm = Elemento('Pm', 'promécio', 3)
Sm = Elemento('Sm', 'samário', 3)
Eu = Elemento('Eu', 'európio', 3)
Gd = Elemento('Gd', 'gadolínio', 3)
Tb = Elemento('Tb', 'térbio', 3)
Dy = Elemento('Dy', 'disprósio', 3)
Ho = Elemento('Ho', 'holmio', 3)
Er = Elemento('Er', 'érbio', 3)
Tm = Elemento('Tm', 'túlio', 3)
Yb = Elemento('Yb', 'itrbio', 3)
Lu = Elemento('Lu', 'lutécio', 3)
Hf = Elemento('Hf', 'háfnio', 4)
Ta = Elemento('Ta', 'tântalo', 5)
W = Elemento('W', 'tungstênio', 6)
Re = Elemento('Re', 'rênio', 7)
Os = Elemento('Os', 'ósmio', 8)
Ir = Elemento('Ir', 'irídio', 4)
Pt = Elemento('Pt', 'platina', 4)
Au = Elemento('Au', 'ouro', 3)
Hg = Elemento('Hg', 'mercúrio', 2)
Tl = Elemento('Tl', 'tálio', 3)
Pb = Elemento('Pb', 'chumbo', 4)
Bi = Elemento('Bi', 'bismuto', 3)
Th = Elemento('Th', 'tório', 4)
Pa = Elemento('Pa', 'protactínio', 5)
U = Elemento('U', 'urânio', 6)
Np = Elemento('Np', 'netúnio', 7)
Pu = Elemento('Pu', 'plutônio', 6)
Am = Elemento('Am', 'amerício', 6)
Cm = Elemento('Cm', 'cúrio', 3)
Bk = Elemento('Bk', 'berquélio', 4)
Cf = Elemento('Cf', 'califórnio', 2)
Es = Elemento('Es', 'einstênio', 3)
Fm = Elemento('Fm', 'férmio', 3)
Md = Elemento('Md', 'mendelévio', 3)
No = Elemento('No', 'nobélio', 3)
Lr = Elemento('Lr', 'laurêncio', 3)
Rf = Elemento('Rf', 'rutherfórdio', 4)
Db = Elemento('Db', 'dúbnio', 5)
Sg = Elemento('Sg', 'seabórgio', 6)
Bh = Elemento('Bh', 'bóhrio', 7)
Hs = Elemento('Hs', 'hássio', 8)
Mt = Elemento('Mt', 'meitnério', 9)
Ds = Elemento('Ds', 'darmstádtio', 10)
Rg = Elemento('Rg', 'roentgênio', 11)
Cn = Elemento('Cn', 'copernício', 12)
Nh = Elemento('Nh', 'nihônio', 13)
Fl = Elemento('Fl', 'fleróvio', 14)
Mc = Elemento('Mc', 'moscóvio', 15)
Lv = Elemento('Lv', 'livermório', 16)
Ts = Elemento('Ts', 'tennessino', 17)
Og = Elemento('Og', 'oganessônio', 18)


def main():
    # Lista de elementos químicos disponíveis
    lista_elementos = [
    H, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar, K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, 
    Ni, Cu, Zn, Ga, Ge, As, Se, Br, Kr, Rb, Sr, Y, Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn, Sb, Te, 
    I, Xe, Cs, Ba, La, Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu, Hf, Ta, W, Re, Os, Ir, Pt,
    Au, Hg, Tl, Pb, Bi, Th, Pa, U, Np, Pu, Am, Cm, Bk, Cf, Es, Fm, Md, No, Lr, Rf, Db, Sg, Bh, Hs, Mt,
    Ds, Rg, Cn, Nh, Fl, Mc, Lv, Ts, Og
    ]

    # Dicionário para mapear os símbolos aos elementos
    dicionario_elementos = {elem.simbolo: elem for elem in lista_elementos}


    while True:
        try:
            system('cls')

            # Pede ao usuário para inserir os símbolos dos elementos
            simbolo1 = input("Digite o símbolo do primeiro elemento: ").strip()
            simbolo2 = input("Digite o símbolo do segundo elemento: ").strip()

            # Encontra os elementos correspondentes aos símbolos inseridos
            try:
                elem1 = dicionario_elementos[simbolo1]
                elem2 = dicionario_elementos[simbolo2]
            except KeyError:
                print("Um ou ambos os símbolos inseridos não correspondem a um elemento conhecido.")
                return

            # Cria a combinação e exibe o nome
            print('\n', Combinacao(elem1, elem2).nome(), end=' ')
            print(f'({Combinacao(elem1, elem2).formula()})')
            input('\nPressione ENTER para continuar . . . ')
        except:
            pass
    

main()