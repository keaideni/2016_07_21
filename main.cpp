#include "Parameter.h"
#include "OP.h"
#include "SuperEnergy.h"
#include "DMRGP.h"
#include "physics.h"



int OP::Max;
int main()
{
	Parameter para;
	para.read();
	
	OP::Max = para.ParticleNo + 1;

	
	//para.D = 200;
	





	//DMRGP DMRG(para);
	//DMRG.fwave.show();

        QWave fwave;

        std::ifstream infile("./Corr/QWave");
        fwave.read(infile);


	std::ofstream Fdata("./result/data", std::ios_base::out | std::ios_base::app);
	CalcuCorr(para.ParticleNo/2+1, fwave, Fdata);
	calcustructure(fwave, para.ParticleNo/2+1, para.ParticleNo, Fdata);
	Fdata.close();


	calcudensity(para.ParticleNo/2+1, fwave, para.ParticleNo);

	system("pause");

	
}