#include "PNL_min.hpp"
#include <gsl/gsl_blas.h>
#include <iostream>
#include <cmath>

using namespace std;


double PNL_funcao::f(const gsl_vector * x){
	return this->funcao(x,this->param);
}
gsl_vector * PNL_funcao::Df(const gsl_vector * x){
	return this->gradiente(x,this->param);
}
double PNL_g(const double t, const gsl_vector * x,
		   const gsl_vector * d, PNL_funcao func){
	double F;
	gsl_vector * y = gsl_vector_alloc(func.dim);
	gsl_blas_dcopy(x,y);
	gsl_blas_daxpy(t,d,y);
	F = func.f(y);
	gsl_vector_free(y);
	return F; // retorna f( x + td )
}

/*! \brief Seção áurea
 *
 *  Método de busca unidimensional utilizando o método
 *  da seção áurea
 *
 * \param x_atual vetor 
 * \param Parameter Parameter description
 * \return Return parameter description
 */

double PNL_busca_unidim_secao_aurea(gsl_vector *x_atual,
							gsl_vector *d, PNL_funcao func,
							double erro, double *param){
	if (x_atual->size != d->size || d->size != func.dim) {
		cout<<"ERRO 0001\n";
		return NAN;
	}

	double theta1 = 1-PHI, theta2 = PHI, a, s = param[0], b=2*s;
	double u,v, F1, F2;
	F1 = PNL_g(b,x_atual,d,func);
	F2 = PNL_g(s,x_atual,d,func);
	
	while ( F1 < F2) {
		a=s; s=b; F2 = F1; b=2*s;
		F1 = PNL_g(b,x_atual,d,func);
	}
	
	u=a+theta1*(b-a); v=a+theta2*(b-a);
	F1 = PNL_g(u,x_atual,d,func);
	F2 = PNL_g(v,x_atual,d,func);
	while ((b-a)>erro) {
		if ( F1 < F2 ) {
			b = v; v = u; u = a + theta1*(b-a);
			F2 = F1; F1 = PNL_g(u,x_atual,d,func);
		} else {
			a = u; u = v; v = a + theta2*(b-a);
			F1 = F2; F2 = PNL_g(v,x_atual,d,func); 
		}
	}

	return (u+v)/2;

}

double PNL_busca_unidim_armijo(gsl_vector *x_atual,
							gsl_vector *d, PNL_funcao func,
							double erro, double *param){
	
	if (x_atual->size != d->size || d->size != func.dim) {
		cout<<"ERRO 0001\n";
		return NAN;
	}
	double F0, F1, Dfd, t=param[0];
	gsl_blas_ddot(func.Df(x_atual),d,&Dfd);
	Dfd *= param[1];
	F0 = func.f(x_atual);
	F1 = PNL_g(t,x_atual,d,func);

	while ( F1 > F0 + t*Dfd) {
		t = param[2]*t;
		F1 = PNL_g(t,x_atual,d,func);
	}
	return t;

}
int PNL_qnewton_rk1(gsl_vector *p, gsl_vector *q, gsl_matrix *H){
	double alpha;
	gsl_vector * v = gsl_vector_alloc(p->size);
	gsl_blas_dcopy(p,v);
	gsl_blas_dsymv(CblasUpper,-1.0,H,q,1.0,v);

	gsl_blas_ddot(v,q,&alpha);
	alpha = 1/alpha;
	gsl_blas_dsyr(CblasUpper,alpha,v,H);
	
	gsl_vector_free(v);
	return 0;
}
int PNL_qnewton_dfp(gsl_vector *p, gsl_vector *q, gsl_matrix *H){
	double alpha, theta;
	gsl_vector * v = gsl_vector_alloc(q->size);

	gsl_blas_ddot(p,q,&alpha);
	alpha = 1/alpha;
	
	gsl_blas_dsymv(CblasUpper,1.0,H,q,0.0,v);

	gsl_blas_ddot(q,v,&theta);
	theta = -1/theta;

	gsl_blas_dsyr(CblasUpper,alpha,p,H);
	gsl_blas_dsyr(CblasUpper,alpha,v,H);

	
	gsl_vector_free(v);
	return 0;
}
int PNL_qnewton_bfgs(gsl_vector *p, gsl_vector *q, gsl_matrix *H){
	double alpha, theta, temp;
	gsl_vector * v = gsl_vector_alloc(q->size);

	gsl_blas_dsymv(CblasUpper,1.0,H,q,0.0,v);

	gsl_blas_ddot(q,q,&temp);
	gsl_blas_ddot(p,q,&alpha);
	theta = -1/temp;
	temp *= alpha;
	gsl_blas_ddot(q,v,&alpha);
	alpha = (1+alpha)/temp;

	
	gsl_blas_dsyr(CblasUpper,alpha,p,H);
	gsl_blas_dger(theta,p,v,H);
	gsl_blas_dger(theta,v,q,H);

	gsl_vector_free(v);
	return 0;
}


void PNL_metodo_desc::iterar(gsl_vector *x_atual, PNL_funcao f, gsl_vector *x_prox, double erro){
	gsl_vector * y = gsl_vector_alloc(f.dim);
	gsl_blas_dsymv(CblasUpper,-1.0,this->H,f.Df(x_atual),0.0,this->d); // d <- -( H * Df )
	this->passo = this->t_passo(x_atual,d,f,erro,this->param);
	cout<<this->passo<<endl;

	gsl_blas_dcopy(x_atual,y);
	gsl_blas_daxpy(this->passo,this->d,y);
	gsl_blas_dcopy(y,x_prox);

	if(this->atualizar_H != NULL){
		gsl_vector *p, *q;
		p = gsl_vector_alloc(f.dim);
		q = gsl_vector_alloc(f.dim);

		gsl_blas_dcopy(x_prox,p);
		gsl_blas_daxpy(-1.0,x_atual,p);

		gsl_blas_dcopy(f.Df(x_prox),q);
		gsl_blas_daxpy(-1.0,f.Df(x_atual),q);

		this->atualizar_H(p,q,this->H);
		gsl_vector_free(p);
		gsl_vector_free(q);
	}
	gsl_vector_free(y);
}


PNL_metodo_desc::PNL_metodo_desc(char tipo, char tbusca, double *param, int dim){
	this->passo = 0;
	this->param = param;
	this->d = gsl_vector_alloc(dim);
	gsl_vector_set_zero(this->d);
	this->H = gsl_matrix_alloc(dim,dim);
	gsl_matrix_set_identity(this->H);
	if (tipo == PNL_QNEWTON_RK1) {
		this->atualizar_H = PNL_qnewton_rk1;
	} else if (tipo == PNL_QNEWTON_DFP) {
		this->atualizar_H = PNL_qnewton_dfp;
	} else if (tipo == PNL_QNEWTON_BFGS) {
		this->atualizar_H = PNL_qnewton_bfgs;
	} else this->atualizar_H = NULL;

	if (tbusca == PNL_UD_SEC_AUREA) {
		this->t_passo = PNL_busca_unidim_secao_aurea;
	} else if (tbusca == PNL_UD_ARMIJO) {
		this->t_passo = PNL_busca_unidim_armijo;
	}
}

void PNL_metodo_desc::def_param(double * param){
}


PNL_min_ambiente::PNL_min_ambiente(char tipo, char t_busca_unidim, int dim){
	double * param;
	this->tipo = tipo;
	this->tbusca = t_busca_unidim;
	gsl_vector *x_ini = gsl_vector_alloc(dim);
	gsl_vector *x_atual = gsl_vector_alloc(dim);
	gsl_vector *x_ant = gsl_vector_alloc(dim);
	gsl_vector *grad = gsl_vector_alloc(dim);

	if (t_busca_unidim == PNL_UD_SEC_AUREA) {
		param = new double(1);
		*param = PNL_SEC_AUREA_RO;
	} else if (t_busca_unidim == PNL_UD_ARMIJO) {
		param = new double(3);
		param[0] = PNL_ARMIJO_TBAR;
		param[1] = PNL_ARMIJO_ETA;
		param[2] = PNL_ARMIJO_THETA;
	}

	this->metodo = new PNL_metodo_desc(tipo,t_busca_unidim, param, dim);
}

PNL_min_ambiente::PNL_min_ambiente(char tipo,
								char t_busca_unidim,
								PNL_funcao f,
								gsl_vector *x_inicial,
								double erro){
	double * param;
	this->tipo = tipo;
	this->tbusca = t_busca_unidim;
	gsl_vector *x_ini = gsl_vector_alloc(f.dim);
	gsl_vector *x_atual = gsl_vector_alloc(f.dim);
	gsl_vector *x_ant = gsl_vector_alloc(f.dim);
	gsl_vector *grad = gsl_vector_alloc(f.dim);

	if (t_busca_unidim == PNL_UD_SEC_AUREA) {
		param = new double(1);
		*param = PNL_SEC_AUREA_RO;
	} else if (t_busca_unidim == PNL_UD_ARMIJO) {
		param = new double(3);
		param[0] = PNL_ARMIJO_TBAR;
		param[1] = PNL_ARMIJO_ETA;
		param[2] = PNL_ARMIJO_THETA;
	}

	this->metodo = new PNL_metodo_desc(tipo,t_busca_unidim, param, f.dim);
	this->def(f,x_inicial,erro);
}

void PNL_min_ambiente::def(PNL_funcao f, gsl_vector *x_inicial, double erro){
	this->funcao = f;
	this->x_ini = x_inicial;
	this->x_ant = gsl_vector_alloc(f.dim);
	this->x_atual = gsl_vector_alloc(f.dim);
	this->grad = gsl_vector_alloc(f.dim);
	gsl_blas_dcopy(this->x_ini,this->x_atual);
	gsl_blas_dcopy(this->x_ini,this->x_ant);
	gsl_blas_dcopy(f.Df(x_inicial),grad);
	this->imag_atual = f.f(x_inicial);
	this->imag_ant = this->imag_atual;
	this->erro = erro;
}

void PNL_min_ambiente::passo(){
	gsl_blas_dcopy(x_atual,x_ant);
	imag_ant = imag_atual;
	this->metodo->iterar(x_ant,funcao,x_atual,erro);
	gsl_blas_dcopy(funcao.Df(x_atual),grad);
	imag_atual = funcao.f(x_atual);
}
int PNL_min_ambiente::teste_gradiente(){
	if (gsl_blas_dnrm2(this->grad) > this->erro) return 1;
	else return 0;
}

void PNL_min_ambiente::calcular(){
	while (teste_gradiente()) {
		this->passo();
	}
}
