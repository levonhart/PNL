#include <PNL_min.hpp>
#include <iostream>

double f_quad (const gsl_vector *x, double * param){
	return param[0]*gsl_vector_get(x, 0)*gsl_vector_get(x, 0) + param[1]*gsl_vector_get(x, 1)*gsl_vector_get(x, 1);
}

gsl_vector * grad_f_quad(const gsl_vector *x, double * param){
	gsl_vector * y = gsl_vector_alloc(2);
	gsl_vector_set(y, 0, 2.0*param[0]*gsl_vector_get(x, 0));
	gsl_vector_set(y, 1, 2.0*param[1]*gsl_vector_get(x, 1));
	return y;
}
void grad_f_quad2(const gsl_vector *x, double *param, gsl_vector * grad){
	grad = grad_f_quad(x,param);
}


int main(int argc, char *argv[])
{
	PNL_funcao func;
	func.param = new double(2);
	func.param[0] = func.param[1] = 4;
	func.funcao = f_quad;
	func.gradiente = grad_f_quad;
	func.dim = 2;

	gsl_vector * x = gsl_vector_alloc(2);
	gsl_vector_set_all(x,1);

	// PNL_min_ambiente * wp = new PNL_min_ambiente(PNL_GRAD,PNL_UD_SEC_AUREA,func,x,0.01);
	PNL_min_ambiente * wp = new PNL_min_ambiente(PNL_QNEWTON_BFGS,PNL_BL_SEC_AUREA,func,x,0.01);
	// PNL_min_ambiente * wp = new PNL_min_ambiente(PNL_GRAD,PNL_UD_SEC_AUREA,2);
	while (wp->teste_gradiente()) {
		printf("x_atual = [ %lf\t%lf ]\tf(x) = %5.10lf\n\n",\
				gsl_vector_get(wp->x_atual,0),gsl_vector_get(wp->x_atual,1),\
				 wp->imag_atual);
		wp->passo();
		
	}
	printf("x_min  = [ %lf\t%lf ]\tf(x) = %5.10lf\n\n",\
			gsl_vector_get(wp->x_atual,0),gsl_vector_get(wp->x_atual,1),\
			 wp->imag_atual);


	return 0;
}
