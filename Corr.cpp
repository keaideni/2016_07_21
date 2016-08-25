#include "Corr.h"

Corr::Corr(){}

Corr::Corr(const Corr& corr)
{
        Orbital = corr.Orbital;

        CorrO = corr.CorrO;

        Type = corr.Type;

}




void Corr::Initial(const Sub& sub,const Sub& subm, const int& orbital, const int& type, const OP& truncU_)
{
        clear();
        Orbital = orbital;
        if(type == 4)
        {
                Type = 3;
        }else
        {
                Type = type;
        }
        OP tempCorrO;
        if(type == 1)
        {
                tempCorrO = subm.SubSysC;

        }else if(type == 2)
        {
                tempCorrO = subm.SubSysCdag;
        }else if(type == 3)//initial for ni which except n1
        {
                tempCorrO.time(subm.SubSysCdag, subm.SubSysC);
        }else if(type == 4)//initial for the n1
        {
                tempCorrO.time(sub.SubSysCdag, sub.SubSysC);
        }else if(type == 5)//initial for the n1n1
        {
                OP tempCorrO1;
                tempCorrO1.time(sub.SubSysCdag, sub.SubSysC);
                tempCorrO.time(tempCorrO1, tempCorrO1);
        }else if(type == 6)//initialize the n1n2
        {
                OP tempCorrO1;
                tempCorrO1.time(sub.SubSysCdag, sub.SubSysC);
                tempCorrO.time(subm.SubSysCdag, subm.SubSysC);
                CorrO.kronO(tempCorrO1, tempCorrO);
                Type = 5;

        }else if(type == 7)//initialize the n2n2
        {
                OP tempCorrO1;
                tempCorrO1.time(subm.SubSysCdag, subm.SubSysC);
                tempCorrO.time(tempCorrO1, tempCorrO1);
        }else if(type == 8)//initialize the n_{l/2+1}
        {
                tempCorrO.time(sub.SubSysCdag, sub.SubSysC);
                Type = 3;
        }else if(type == 9)//initialize the Cdag_{l/2+1}.l=para.latticesize
        {
                tempCorrO=sub.SubSysC;
                Type = 1;
        }else if(type == 10)
        {
                tempCorrO=sub.SubSysCdag;
                Type = 2;
        }
        if(type == 4 || type == 5 || type == 8 || type==9||type==10)
        {
                CorrO.kronO(tempCorrO, subm.SubSysEye);
        }else if(type == 6)
        {

        }else
        {
                CorrO.kronO(sub.SubSysEye, tempCorrO);
        }
        CorrO.trunc(truncU_);


}

//initialize the n1ni and n2ni.
void Corr::Initial(const Corr& corr, const Sub& m, const int& orbital, const int& type, const OP& truncU)
{

        clear();
        Orbital = orbital;
        Type = type;

        OP tempCorrO;
        tempCorrO.time(m.SubSysCdag, m.SubSysC);
        CorrO.kronO(corr.CorrO, tempCorrO);

        CorrO.trunc(truncU);
}




OP Corr::corro()
{
        return CorrO;
}




int Corr::orbital()
{
        return Orbital;
}



int Corr::type()
{
        return Type;
}



void Corr::update(const Sub& m, const OP& truncU, const int& way)
{
        if(way == 1)
        {
               OP tempCorrO;
               tempCorrO.kronO(this->CorrO, m.SubSysEye);
               tempCorrO.trunc(truncU); 
               CorrO = tempCorrO;
        }else
        {
                OP tempCorrO;
                tempCorrO.kronO(m.SubSysEye, this->CorrO);

                tempCorrO.trunc(truncU); 
                CorrO = tempCorrO;
        }
}





void Corr::save()
{
        std::string str = itos(Orbital);
        if(Type == 1)
        {
                str = "./Corr/C" + str;
        }else if(Type == 2)
        {
                str = "./Corr/CDag" + str;
        }else if(Type == 3)//save the ni
        {
                str = "./Corr/N" + str;
        }else if(Type == 5)
        {
                str = "./Corr/N1N" + str;
        }else
        {
                str = "./Corr/N2N" + str;
        }

        std::ofstream outfile(str);
        
        outfile << Orbital <<std::endl;
        outfile << Type <<std::endl;
        CorrO.save(outfile);

        outfile.close();
}



void Corr::read(const int& orbital, const int& type)
{
        clear();
        std::string str = itos(orbital);
        if(type == 1)
        {
                str = "./Corr/C" + str;
        }else if(type == 2)
        {
                str = "./Corr/CDag" + str;
        }else if(type == 3)
        {
                str = "./Corr/N" + str;
        }else if(type == 5)
        {
                str = "./Corr/N1N" + str;
        }else
        {
                str = "./Corr/N2N" + str;
        }

        std::ifstream infile(str);
        if(!infile.is_open())
        {
                std::cout<<"the file "<<str<<" doesn't exist!"<<std::endl;
                exit(1);
        }
        infile >> Orbital;
        infile >> Type;

        CorrO.read(infile);


}




void Corr::show()
{
        std::string str = itos(Orbital);
        if(Type == 1)
        {
                str = "the Operator C on the msite "+str;
        }else if(Type == 2)
        {
                str = "the Operator CDag on the msite "+str;
        }else
        {
                str = "the Operator N on the msite "+str;
        }
        std::cout << str <<std::endl;
        CorrO.show();

}



void Corr::clear()
{
        CorrO.clear();
}




