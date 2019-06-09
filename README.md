# PNL
minimizador multidimensional

# Requisitos
  - GSL (GNU Scientific Library)
  - pkg-config (opcional)

# Instruções de uso
Essa biblioteca consiste em dois arquivos, `PNL_min.hpp` e `libpnl.a`. Para utilizar, adicione o seguinte código no 
cabeçalho do seu programa:
```
#include <PNL_min.hpp>
```

Para compilar, utilize a tag `-pnl` e o `pkg-config` com os parametros do GSL.
Ex.:
```
g++ main.cpp -o saida.o -I. -L. -lpnl `pkg-config --libs gsl`
```
Se não possuir a ferramenta `pkg-config`:
```
g++ main.cpp -o saida.o -I. -L. -lpnl -lgsl -lgslcblas -lm
```

Os métodos são realizados no ambiente de minimização `PNL_min_anbiente` sobre a classe `PNL_funcao`.

Para encontrar um mínimo local em uma função de n variáveis definida como
`double minha_função( gsl_vector * x, double * param );`
e com gradiente definido como
`gsl_vector * minha_função_grad( gsl_vector * x, double * param );`
defina:
```
int n = 2 // função de 2 variáveis
PNL_funcao func;
func.param = new double(2);
func.param[0] = func.param[1] = 4;
func.funcao = minha_funcao;
func.gradiente = minha_funcao_grad;
func.dim = n;
```
defina ainda o ponto inicial de busca:
```
gsl_vector * x = gsl_vector_alloc(n);
gsl_vector_set(x,0,2);
gsl_vector_set(x,1,5); // nesse exemplo x = [ 2   5 ]
```

em seguida, inicie o ambiente de minimização
```
PNL_min_ambiente * amb = new PNL_min_ambiente( metodo , metodo_busca_linear , func, x, Erro);
```
aqui, `Erro` é qualquer constante positiva não nula.
`metodo` pode ser:
+ PNL_GRAD - método de descida utilizando o vetor gradiente como direção;
+ PNL_QNEWTON_RK1 - método quase-Newton Rank One, que utiliza a fórmula de rank one correction;
+ PNL_QNEWTON_DFP - método quase-Newton de Davidon-Fletcher-Powell;
+ PNL_QNEWTON_BFGS - método quase-Newton de Broyden-Fletcher-Goldfard-Shanno.

e `metodo_busca_linear pode ser:
+ PNL_BL_SEC_AUREA - método de busca linear utilizando a regra da seção áurea;
+ PNL_BL_ARMIJO - método de busca linear utilizando o critério de Armijo.


por fim, basta utilizar
`amb->calcular();` 
ou, se houver interesse em cada passo da busca, o código abaixo imprime os valores de x e suas respectivas imagem a cada iteração;
```
while (amb->teste_gradiente()) { //testa se o gradinte é nulo, tem norma menor que Erro;
  printf("x_atual = [ %lf\t%lf ]\tf(x) = %5.10lf\n\n",\
			gsl_vector_get(amb->x_atual,0),gsl_vector_get(amb->x_atual,1),\
			 amb->imag_atual);
       
	amb->passo();
}
printf("x_min = [ %lf\t%lf ]\tf(x) = %5.10lf\n\n",\
		gsl_vector_get(amb->x_atual,0),gsl_vector_get(amb->x_atual,1),\
		amb->imag_atual);
```
