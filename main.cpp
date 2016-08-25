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
	





	DMRGP DMRG(para);
	//DMRG.fwave.show();

        QWave fwave;

        std::ifstream infile("./Corr/QWave");
        fwave.read(infile);


	std::ofstream Fdata("./result/data", std::ios_base::out | std::ios_base::app);
	CalcuCorr(DMRG.OrbitalM, fwave, Fdata);
	calcustructure(fwave, DMRG.OrbitalM, para.LatticeSize/2, Fdata);
	Fdata.close();


	calcudensity(DMRG.OrbitalM, fwave, para.LatticeSize/2);

	system("pause");

	
}