#include "stdafx.h"
#include "GeneticAlgorithm.h"

void GeneticAlgorithm::isotropicParameters(double * pOld, double * pNew, double * rand1, double * rand2, double * rand3){
	
		double p = randomDouble(0,1);
		if(p > _crossingProbability){
			pNew[0] = pOld[0];
		}
		else{
			pNew[0] = rand1[0] + _scaleFactor*(rand2[0] - rand3[0]);
		}
		p = randomDouble(0,1);
		if(p > _crossingProbability){
			pNew[3] = pOld[3];
		}
		else{
			pNew[3] = rand1[3] + _scaleFactor*(rand2[3] - rand3[3]);
		}
			pNew[1] = pNew[0];
			pNew[2] = pNew[0];
			pNew[4] = pNew[3];
			pNew[5] = pNew[3];
			pNew[6] = pNew[0] - 2 * pNew[3];
			pNew[7] = pNew[0] - 2 * pNew[3];
			pNew[8] = pNew[0] - 2 * pNew[3];

		return;
}

void GeneticAlgorithm::cubicParameters(double * pOld, double * pNew, double * rand1, double * rand2, double * rand3){
	
		double p = randomDouble(0,1);
		if(p > _crossingProbability){
			pNew[0] = pOld[0];
		}
		else{
			pNew[0] = rand1[0] + _scaleFactor*(rand2[0] - rand3[0]);
		}
		p = randomDouble(0,1);
		if(p > _crossingProbability){
			pNew[3] = pOld[3];
		}
		else{
			pNew[3] = rand1[3] + _scaleFactor*(rand2[3] - rand3[3]);
		}
		p = randomDouble(0,1);
		if(p > _crossingProbability){
			pNew[6] = pOld[6];
		}
		else{
			pNew[6] = rand1[6] + _scaleFactor*(rand2[6] - rand3[6]);
		}
			pNew[1] = pNew[0];
			pNew[2] = pNew[0];
			pNew[4] = pNew[3];
			pNew[5] = pNew[3];
			pNew[7] = pNew[6];
			pNew[8] = pNew[6];

		return;
}

void GeneticAlgorithm::hexagonalParameters(double * pOld, double * pNew, double * rand1, double * rand2, double * rand3){
	
		double p = randomDouble(0,1);
		if(p > _crossingProbability){
			pNew[0] = pOld[0];
		}
		else{
			pNew[0] = rand1[0] + _scaleFactor*(rand2[0] - rand3[0]);
		}
		p = randomDouble(0,1);
		if(p > _crossingProbability){
			pNew[2] = pOld[2];
		}
		else{
			pNew[2] = rand1[2] + _scaleFactor*(rand2[2] - rand3[2]);
		}
		p = randomDouble(0,1);
		if(p > _crossingProbability){
			pNew[3] = pOld[3];
		}
		else{
			pNew[3] = rand1[3] + _scaleFactor*(rand2[3] - rand3[3]);
		}
			p = randomDouble(0,1);
		if(p > _crossingProbability){
			pNew[6] = pOld[6];
		}
		else{
			pNew[6] = rand1[6] + _scaleFactor*(rand2[6] - rand3[6]);
		}
			p = randomDouble(0,1);
		if(p > _crossingProbability){
			pNew[7] = pOld[7];
		}
		else{
			pNew[7] = rand1[7] + _scaleFactor*(rand2[7] - rand3[7]);
		}
			pNew[1] = pNew[0];
			pNew[4] = pNew[3];
			pNew[8] = pNew[7];
			pNew[5] = 0.5*(pNew[0] - pNew[6]);
			
		return;
}

void GeneticAlgorithm::tetragonalParameters(double * pOld, double * pNew, double * rand1, double * rand2, double * rand3){
	
		double p = randomDouble(0,1);
		if(p > _crossingProbability){
			pNew[0] = pOld[0];
		}
		else{
			pNew[0] = rand1[0] + _scaleFactor*(rand2[0] - rand3[0]);
		}

		p = randomDouble(0,1);
		if(p > _crossingProbability){
			pNew[2] = pOld[2];
		}
		else{
			pNew[2] = rand1[2] + _scaleFactor*(rand2[2] - rand3[2]);
		}

		p = randomDouble(0,1);
		if(p > _crossingProbability){
			pNew[3] = pOld[3];
		}
		else{
			pNew[3] = rand1[3] + _scaleFactor*(rand2[3] - rand3[3]);
		}

		p = randomDouble(0,1);
		if(p > _crossingProbability){
			pNew[5] = pOld[5];
		}
		else{
			pNew[5] = rand1[5] + _scaleFactor*(rand2[5] - rand3[5]);
		}

			p = randomDouble(0,1);
		if(p > _crossingProbability){
			pNew[6] = pOld[6];
		}
		else{
			pNew[6] = rand1[6] + _scaleFactor*(rand2[6] - rand3[6]);
		}

			p = randomDouble(0,1);
		if(p > _crossingProbability){
			pNew[7] = pOld[7];
		}
		else{
			pNew[7] = rand1[7] + _scaleFactor*(rand2[7] - rand3[7]);
		}

			pNew[1] = pNew[0];
			pNew[4] = pNew[3];
			pNew[8] = pNew[7];
		
		return;
}

void GeneticAlgorithm::orthorhombicParameters(double * pOld, double * pNew, double * rand1, double * rand2, double * rand3){
	
		double p = randomDouble(0,1);
		if(p > _crossingProbability){
			pNew[0] = pOld[0];
		}
		else{
			pNew[0] = rand1[0] + _scaleFactor*(rand2[0] - rand3[0]);
		}

		p = randomDouble(0,1);
		if(p > _crossingProbability){
			pNew[1] = pOld[1];
		}
		else{
			pNew[1] = rand1[1] + _scaleFactor*(rand2[1] - rand3[1]);
		}

		p = randomDouble(0,1);
		if(p > _crossingProbability){
			pNew[2] = pOld[2];
		}
		else{
			pNew[2] = rand1[2] + _scaleFactor*(rand2[2] - rand3[2]);
		}

		p = randomDouble(0,1);
		if(p > _crossingProbability){
			pNew[3] = pOld[3];
		}
		else{
			pNew[3] = rand1[3] + _scaleFactor*(rand2[3] - rand3[3]);
		}

		p = randomDouble(0,1);
		if(p > _crossingProbability){
			pNew[4] = pOld[4];
		}
		else{
			pNew[4] = rand1[4] + _scaleFactor*(rand2[4] - rand3[4]);
		}

		p = randomDouble(0,1);
		if(p > _crossingProbability){
			pNew[5] = pOld[5];
		}
		else{
			pNew[5] = rand1[5] + _scaleFactor*(rand2[5] - rand3[5]);
		}

			p = randomDouble(0,1);
		if(p > _crossingProbability){
			pNew[6] = pOld[6];
		}
		else{
			pNew[6] = rand1[6] + _scaleFactor*(rand2[6] - rand3[6]);
		}

			p = randomDouble(0,1);
		if(p > _crossingProbability){
			pNew[7] = pOld[7];
		}
		else{
			pNew[7] = rand1[7] + _scaleFactor*(rand2[7] - rand3[7]);
		}

				p = randomDouble(0,1);
		if(p > _crossingProbability){
			pNew[8] = pOld[8];
		}
		else{
			pNew[8] = rand1[8] + _scaleFactor*(rand2[8] - rand3[8]);
		}

		
		return;
}