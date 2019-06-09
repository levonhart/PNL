#ifndef PNL_MIN_HPP_ATD1WBYL
#define PNL_MIN_HPP_ATD1WBYL


#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#define PNL_GRAD 0x47
#define PNL_QNEWTON_RK1 0x52
#define PNL_QNEWTON_DFP 0x44
#define PNL_QNEWTON_BFGS 0x42
#define PNL_BL_ARMIJO 0x41
#define PNL_BL_SEC_AUREA 0x61

#define PHI 0.618033988749894848204586834365638117720309179805762862135448622
#define PNL_SEC_AUREA_RO 1.0
#define PNL_ARMIJO_TBAR 5.0
#define PNL_ARMIJO_ETA 0.5
#define PNL_ARMIJO_THETA 0.5


struct PNL_funcao {
	size_t dim;
	double * param;
	double (*funcao) (const gsl_vector * x, double * param);
	gsl_vector * (*gradiente) (const gsl_vector * x, double * param);
	double f (const gsl_vector * x);
	gsl_vector * Df (const gsl_vector * x);
	void avalie(gsl_vector * x, double * imagem, gsl_vector * gradiente);
};

double PNL_g(const double t, const gsl_vector * x,
		   const gsl_vector * d, PNL_funcao func);

double PNL_busca_unidim_secao_aurea(gsl_vector *x_atual,
							gsl_vector *d, PNL_funcao func,
							double erro, double * param);
double PNL_busca_unidim_secao_aurea(gsl_vector *x_atual,
							gsl_vector *d, PNL_funcao func,
							double erro, double * param);
int PNL_qnewton_rk1(gsl_vector *p, gsl_vector *q, gsl_matrix *H);
int PNL_qnewton_dfp(gsl_vector *p, gsl_vector *q, gsl_matrix *H);
int PNL_qnewton_bfgs(gsl_vector *p, gsl_vector *q, gsl_matrix *H);

struct PNL_metodo_desc {
	double passo, *param;
	gsl_matrix *H;
	gsl_vector *d;
	double (*t_passo) (gsl_vector *x, gsl_vector *dir,
			PNL_funcao f, double erro, double *param);
	int (*atualizar_H)(gsl_vector * p, gsl_vector * q, gsl_matrix * H);

	void iterar (gsl_vector *x_atual, PNL_funcao f, gsl_vector *x_prox, double erro);
	PNL_metodo_desc(char tipo, char tbusca, double * param, int dim);
	void def_param(double * param);

};

struct PNL_min_ambiente {
	char tipo, tbusca;
	gsl_vector *x_ini, *x_atual, *x_ant, *grad;
	double imag_atual, imag_ant, erro;
	PNL_funcao funcao;
	PNL_metodo_desc * metodo;

	PNL_min_ambiente(char tipo, char t_busca_linear, int dim);
	PNL_min_ambiente(char tipo,
					char t_busca_linear,
					PNL_funcao f,
					gsl_vector *x_inicial,
					double erro);

	void def(PNL_funcao f, gsl_vector *x_inicial, double erro);
	void passo();
	void calcular();
	int teste_gradiente();

};



#endif /* end of include guard: PNL_MIN_HPP_ATD1WBYL */
