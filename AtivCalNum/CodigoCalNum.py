# Codigo que que realzia as tarefas da atividade de Calculo Numérico
#Bibliotecas importadas:
import numpy as np

#---------------------------------Funções da Interpolação Polinomial---------------------------------
# Função de Interpolação Polinomial:
def interpolacao_polinomial(x, y):
    # Calculando os coeficientes do polinômio interpolador:
    coef_poli_inter = np.polyfit(x, y, 1)

    # Tratando erros, arredondando o coeficiente do polinômio interpolador:
    coef_poli_inter = np.round(coef_poli_inter, 5)

    # Criando o Polinômio Interpolador:
    #polinomio = np.poly1d(coef_poli_inter)

    return(coef_poli_inter)

# Função para aplicar o Polinômio Interpolador em novas entradas
def aplicar_inter_poli(array_poli, new_values):
    # Criando um array auxiliar
    aux_resultado = np.array([])

    for number in new_values:
        i = 0
        aux_poli = []
        while(i < len(array_poli)):
            aux_poli.append(number**(len(array_poli) - 1 - i) * (array_poli[i]))
            i = i + 1
        aux_resultado = np.append(aux_resultado, sum(aux_poli))
    return aux_resultado
            

#---------------------------------Funções para Calcular a Integral---------------------------------
# Função para Calcular a área do Trapezio:
def areaTrapezio(b, B, h):
    return ((B + b)*h/2)

# Função para gerar a lista de Áreas dos Trapézios
def gerarListaAreas(listX, listY):
    listAreas = [];
    i = 0;
    while(i < len(listX) - 1):
        listAreas.append(areaTrapezio(listY[i+1], listY[i], listX[i+1] - listX[i]))
        i += 1

    return listAreas

#---------------------------------Funções do Ajuste de Curvas---------------------------------
# Função que calcula o Coeficiente de Inclinação:
def CoefInclinacao(x, y, x2, xy):
    return (len(x) * sum(xy) - sum(x) * sum(y)) / (len(x) * sum(x2) - sum(x) * sum(x))
    
# Função que calcula o Intercepto:
def Intercepto(x, y, CoefInc):
    return  (sum(y) / len(y)) - (sum(x) / len(x)) * CoefInc
#---------------------------------Funções do Ajuste de Curvas---------------------------------
# Função f(x) = 2x
def funcX(x):
    return (x*x*x - x - 2)

# Função que realiza o metodo da Bisseção 
def bissecao(a, b, tol):
    # Checa se é possível realizar o método:
    if funcX(a) * funcX(b) >= 0:
        print(f'O método da Bisseção não é aplicável!')
        return None
        
    while (b-a)/2 > tol:
        # Calculo do ponto médio:
        c = (a + b)/2 
        if(funcX(c) == 0):
            # Encontrou a raiz exata, logo retorna o c!
            return c
        elif(funcX(a) * funcX(c) < 0):
            # Atualiza o intervalo para [b,c]!
            b = c
        else:
            # Atualiza o intervalo para [c, b]!
            a = c
    return (a + b)/2
#--------------------------------------MAIN--------------------------------------

# Declaração do X:
x = np.array([0,2,4,6,8,10])

 # Declaração do Y:
y = np.array([0,4,8,12,16,20])

# ---------- Interpolação Polinomial ----------
tempos_n_medidos = np.array([1,3,5,7,9,11])

print(f"Interpolação Polinolial: {interpolacao_polinomial(x, y)}")

print(f"Interpolação Polinolial Aplicada: {aplicar_inter_poli(interpolacao_polinomial(x, y), tempos_n_medidos)}")

# ---------- Integral ----------
print(f'\n')
lisAreas = gerarListaAreas(x, y)

# A integral é a soma da lista de Trapézios
print(f'Integral : {sum(lisAreas)}')

#----------Ajuste de Curvas----------

print(f'\n')
# Declarando X * Y e X * X:
XxY = []
XxX = []
i = 0
# Populando X * Y e X * X:
while(i<len(x)):
    XxY.append(x[i]*y[i])
    XxX.append(x[i]*x[i])
    i = i + 1

# Chamando a função que calcula o Coeficiente de Inclinação:
CoefInc = CoefInclinacao(x,y,XxX,XxY)

# Print da função que calcula o Coeficiente de Inclinação:
print(f'Coeficiente de Inclinação: {CoefInc}')

# Print e chamando a função que calcula o Intercepto:
print(f'Intercepto: {Intercepto(x, y, CoefInc)}')

# Print da função formada pelo Interceptor e Coeficiente de Inclinação:
print(f'A função é: f(x) = {Intercepto(x, y, CoefInc)} + ({CoefInc} * x)')

#----------Método da Bisseção----------
print(f'\n')
# Declarando os valores para realizar p Método da Bisseção
Xmin = 0
Xmax = 200
tol = 1e-5

raiz = bissecao(Xmin, Xmax, tol)
# Checa se foi realizado o Método da Bisseção:
if raiz != None:
    print(f'A raiz é {raiz}')